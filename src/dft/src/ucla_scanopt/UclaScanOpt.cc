// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "UclaScanOpt.hh"

#include <atomic>
#include <array>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <numeric>

#include "ortools/graph/assignment.h"

#include <ABKCommon/abkseed.h>
#include <ScanOpt/optimizer.h>
#include <ScanOpt/scanTourDZ.h>

namespace dft {
namespace {

class FlatScanChain : public abkscanopt::ScanChain
{
 public:
  struct CellData
  {
    int in_x = 0;
    int in_y = 0;
    int out_x = 0;
    int out_y = 0;
    std::string name;
  };

  FlatScanChain(const std::vector<CellData>& cells,
                unsigned begin_idx,
                unsigned end_idx)
      : abkscanopt::ScanChain()
  {
    for (const auto& c : cells) {
      _addCell(c.in_x, c.in_y, c.out_x, c.out_y, c.name);
    }

    _normalizePO();

    // Fix endpoints via legalins/legalouts; Optimizer1 keeps path endpoints
    // fixed to the initial path endpoints.
    _addLegalin(static_cast<int>(begin_idx));
    _addLegalout(static_cast<int>(end_idx));

    _checkIndexRanges();
    _checkInOutLegality();

    _path.clear();
    _path.reserve(cells.size());
    for (unsigned i = 0; i < cells.size(); ++i) {
      _path.push_back(i);
    }
    _pathValid = true;
    computeInverse();
  }

  void restorePath(const std::vector<unsigned>& path)
  {
    _path = path;
    _pathValid = true;
    computeInverse();
  }

  std::vector<unsigned>& mutablePath() { return _path; }
};

unsigned clampToUnsigned(uint64_t v)
{
  return static_cast<unsigned>(
      std::min<uint64_t>(v, std::numeric_limits<unsigned>::max()));
}

unsigned clampToDeterministicSeed(uint64_t v)
{
  const unsigned s = clampToUnsigned(v);
  // ABKCommon uses UINT_MAX as the sentinel for "no explicit seed".
  if (s == std::numeric_limits<unsigned>::max()) {
    return std::numeric_limits<unsigned>::max() - 1;
  }
  return s;
}

void configureSeedHandlerOnce(uint64_t seed)
{
  static std::once_flag seed_once;
  std::call_once(seed_once, [seed]() {
    // UCLApack's ABKCommon SeedHandler writes a `seeds.out` lock/log file in the
    // CWD by default. Disable that so embedded use (and parallel runs) don't
    // create/contend on a global file.
    SeedHandler::turnOffLogging();

    // Force a deterministic external seed so any ABKCommon RNGs created without
    // an explicit seed (or using multipartite locIdent seeds) remain
    // deterministic across runs.
    SeedHandler::overrideExternalSeed(clampToDeterministicSeed(seed));
  });
}

double optimizeLevelSafe(FlatScanChain& chain,
                         RandomRawUnsigned& randuns,
                         const abkscanopt::Optimizer::Params& params)
{
  abkscanopt::ScanTourDZ::OptParams opt_params;
  opt_params.nDescents = params.nDescents;
  opt_params.kickMove = params.kickMove;
  opt_params.only2Opt = params.only2Opt;
  opt_params.temp_control = params.temp_control;
  opt_params.bZeroFix = true;

  const unsigned collapse_size = chain.getNumCells() - 1;
  const unsigned in_idx = chain.getPath()[0];
  const unsigned out_idx = chain.getPath()[collapse_size];

  unsigned nnear = params.nnear;
  if (nnear > collapse_size - 1) {
    nnear = collapse_size - 1;
  }

  abkscanopt::ScanTourDZ tour(chain, nnear, randuns, opt_params);
  if (!chain.isSubpath()) {
    tour.scoOptAll();

    unsigned lesser = in_idx;
    unsigned greater = out_idx;
    if (lesser > greater) {
      std::swap(lesser, greater);
    }

    std::vector<unsigned> collapsed_path = tour.getPath();
    auto zero_it = std::find(collapsed_path.begin(), collapsed_path.end(), 0U);
    if (zero_it == collapsed_path.end()) {
      throw std::runtime_error("UclaScanOptOrder: missing distinguished zero");
    }

    // UCLApack's Optimizer1 expects the distinguished-zero cell to remain in
    // position 0 of the tour. Some cases (notably when there are no partial
    // order constraints) can return a rotated tour where the zero cell is not
    // first, which later triggers an ABKCommon fatal error while "uncollapsing"
    // the endpoints. Canonicalize by rotating the cycle so the zero cell is
    // always first.
    if (zero_it != collapsed_path.begin()) {
      std::rotate(collapsed_path.begin(), zero_it, collapsed_path.end());
    }

    std::vector<unsigned>& path = chain.mutablePath();
    for (unsigned i = 1; i < collapse_size; ++i) {
      const unsigned idx = collapsed_path[i];
      if (idx == 0U) {
        throw std::runtime_error(
            "UclaScanOptOrder: internal distinguished zero after rotation");
      }
      path[i] = abkscanopt::ScanCells::uncollapseIndex(idx, lesser, greater);
    }

    chain.computeInverse();
  }

  return tour.scoTourCost();
}

void optimizeChainSafe(FlatScanChain& chain,
                       RandomRawUnsigned& randuns,
                       const abkscanopt::Optimizer::Params& params,
                       const std::optional<std::chrono::steady_clock::time_point>&
                           deadline)
{
  const unsigned major_loops = std::max(1U, params.majorLoops);

  std::vector<unsigned> best_path = chain.getPath();
  double best_cost = std::numeric_limits<double>::infinity();

  std::vector<unsigned> last_path;
  double last_cost = std::numeric_limits<double>::infinity();
  bool have_last = false;

  for (unsigned loop = 0; loop < major_loops; ++loop) {
    if (deadline.has_value()
        && std::chrono::steady_clock::now() >= deadline.value()) {
      break;
    }
    const double cost = optimizeLevelSafe(chain, randuns, params);

    if (cost < best_cost) {
      best_cost = cost;
      best_path = chain.getPath();
    }

    if (params.zeroTemp) {
      if (!have_last || cost < last_cost) {
        have_last = true;
        last_cost = cost;
        last_path = chain.getPath();
      } else {
        chain.restorePath(last_path);
      }
    }
  }

  if (!best_path.empty()) {
    chain.restorePath(best_path);
  }
}

uint64_t splitmix64(uint64_t x)
{
  // http://xorshift.di.unimi.it/splitmix64.c
  uint64_t z = x + 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

void doubleBridgeKick(std::vector<unsigned>& path, std::mt19937_64& rng)
{
  if (path.size() < 8) {
    return;
  }
  const std::size_t n = path.size();
  if (n <= 3) {
    return;
  }

  // Keep begin/end fixed.
  std::uniform_int_distribution<std::size_t> dist(1, n - 2);

  std::array<std::size_t, 4> cuts;
  for (int tries = 0; tries < 20; ++tries) {
    for (std::size_t& c : cuts) {
      c = dist(rng);
    }
    std::sort(cuts.begin(), cuts.end());
    if (cuts[0] < cuts[1] && cuts[1] < cuts[2] && cuts[2] < cuts[3]) {
      break;
    }
  }

  if (!(cuts[0] < cuts[1] && cuts[1] < cuts[2] && cuts[2] < cuts[3])) {
    return;
  }

  const std::size_t a = cuts[0];
  const std::size_t b = cuts[1];
  const std::size_t c = cuts[2];
  const std::size_t d = cuts[3];

  std::vector<unsigned> next;
  next.reserve(n);
  next.insert(next.end(), path.begin(), path.begin() + a);
  next.insert(next.end(), path.begin() + c, path.begin() + d);
  next.insert(next.end(), path.begin() + b, path.begin() + c);
  next.insert(next.end(), path.begin() + a, path.begin() + b);
  next.insert(next.end(), path.begin() + d, path.end());
  path.swap(next);
}

void reverseSegmentKick(std::vector<unsigned>& path, std::mt19937_64& rng)
{
  if (path.size() < 6) {
    return;
  }
  const std::size_t n = path.size();
  // Keep begin/end fixed.
  std::uniform_int_distribution<std::size_t> dist(1, n - 2);
  std::size_t a = dist(rng);
  std::size_t b = dist(rng);
  if (a > b) {
    std::swap(a, b);
  }
  if (b - a < 2) {
    return;
  }
  std::reverse(path.begin() + a, path.begin() + b);
}

bool orderCrossoverOX(const std::vector<unsigned>& parent_a,
                      const std::vector<unsigned>& parent_b,
                      std::mt19937_64& rng,
                      std::vector<unsigned>& out_path)
{
  if (parent_a.size() != parent_b.size()) {
    return false;
  }
  const std::size_t n = parent_a.size();
  if (n < 6) {
    return false;
  }

  // Endpoints must match so we can keep them fixed.
  if (parent_a.front() != parent_b.front() || parent_a.back() != parent_b.back()) {
    return false;
  }

  // Choose an interior segment [a, b] (positions, not node IDs).
  std::uniform_int_distribution<std::size_t> dist(1, n - 2);
  std::size_t a = dist(rng);
  std::size_t b = dist(rng);
  if (a > b) {
    std::swap(a, b);
  }
  if (a == b) {
    return false;
  }

  // Avoid extremely short segments; fall back to mutation if we can't pick.
  if (b - a < 2) {
    return false;
  }

  const unsigned kUnset = std::numeric_limits<unsigned>::max();
  out_path.assign(n, kUnset);
  out_path.front() = parent_a.front();
  out_path.back() = parent_a.back();

  std::vector<char> used(n, 0);
  used[out_path.front()] = 1;
  used[out_path.back()] = 1;

  for (std::size_t pos = a; pos <= b; ++pos) {
    const unsigned node = parent_a[pos];
    if (node >= n) {
      return false;
    }
    out_path[pos] = node;
    used[node] = 1;
  }

  // Fill remaining positions in cyclic order starting at b+1.
  std::size_t fill_pos = b + 1;
  if (fill_pos >= n - 1) {
    fill_pos = 1;
  }
  auto advance = [&]() {
    ++fill_pos;
    if (fill_pos >= n - 1) {
      fill_pos = 1;
    }
  };

  for (std::size_t pos = 1; pos + 1 < n; ++pos) {
    const unsigned node = parent_b[pos];
    if (node >= n || used[node]) {
      continue;
    }
    while (out_path[fill_pos] != kUnset) {
      advance();
    }
    out_path[fill_pos] = node;
    used[node] = 1;
    advance();
  }

  for (std::size_t pos = 1; pos + 1 < n; ++pos) {
    if (out_path[pos] == kUnset) {
      return false;
    }
  }

  return true;
}

int64_t uclaEdgeCost(const FlatScanChain::CellData& a,
                     const FlatScanChain::CellData& b)
{
  return static_cast<int64_t>(std::abs(a.out_x - b.in_x))
         + static_cast<int64_t>(std::abs(a.out_y - b.in_y));
}

int64_t uclaPathCost(const std::vector<FlatScanChain::CellData>& cells,
                     const std::vector<unsigned>& path)
{
  int64_t cost = 0;
  for (std::size_t i = 0; i + 1 < path.size(); ++i) {
    const auto& a = cells[path[i]];
    const auto& b = cells[path[i + 1]];
    cost += uclaEdgeCost(a, b);
  }
  return cost;
}

struct UclaPathScore
{
  int64_t cost = std::numeric_limits<int64_t>::max();
  int64_t max_edge = std::numeric_limits<int64_t>::max();
  int64_t tail_sum = std::numeric_limits<int64_t>::max();
};

UclaPathScore uclaEvaluatePath(const std::vector<FlatScanChain::CellData>& cells,
                               const std::vector<unsigned>& path)
{
  UclaPathScore score;
  if (path.size() < 2) {
    score.cost = 0;
    score.max_edge = 0;
    score.tail_sum = 0;
    return score;
  }

  constexpr std::size_t kTailEdges = 4;
  std::array<int64_t, kTailEdges> worst_edges{};
  score.cost = 0;
  score.max_edge = 0;
  score.tail_sum = 0;

  for (std::size_t i = 0; i + 1 < path.size(); ++i) {
    const auto& a = cells[path[i]];
    const auto& b = cells[path[i + 1]];
    const int64_t edge = uclaEdgeCost(a, b);
    score.cost += edge;
    score.max_edge = std::max(score.max_edge, edge);

    for (std::size_t pos = 0; pos < worst_edges.size(); ++pos) {
      if (edge <= worst_edges[pos]) {
        continue;
      }
      for (std::size_t shift = worst_edges.size(); shift > pos + 1; --shift) {
        worst_edges[shift - 1] = worst_edges[shift - 2];
      }
      worst_edges[pos] = edge;
      break;
    }
  }

  for (const int64_t edge : worst_edges) {
    score.tail_sum += edge;
  }
  return score;
}

bool uclaScoreLess(const UclaPathScore& a,
                   const std::vector<unsigned>& a_path,
                   const UclaPathScore& b,
                   const std::vector<unsigned>& b_path)
{
  if (a.cost != b.cost) {
    return a.cost < b.cost;
  }
  if (a.max_edge != b.max_edge) {
    return a.max_edge < b.max_edge;
  }
  if (a.tail_sum != b.tail_sum) {
    return a.tail_sum < b.tail_sum;
  }
  return std::lexicographical_compare(
      a_path.begin(), a_path.end(), b_path.begin(), b_path.end());
}

bool uclaMetricLess(const UclaPathScore& a, const UclaPathScore& b)
{
  if (a.cost != b.cost) {
    return a.cost < b.cost;
  }
  if (a.max_edge != b.max_edge) {
    return a.max_edge < b.max_edge;
  }
  return a.tail_sum < b.tail_sum;
}

bool uclaTimeExpired(
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  return deadline.has_value()
         && std::chrono::steady_clock::now() >= deadline.value();
}

bool uclaAdjacentSwapImprove(const std::vector<FlatScanChain::CellData>& cells,
                             std::vector<unsigned>& order,
                             bool end_fixed)
{
  bool improved = false;
  const std::size_t n = order.size();
  if (n < 3) {
    return false;
  }
  const std::size_t limit = end_fixed ? (n - 1) : n;
  for (std::size_t pos = 1; pos + 1 < limit; ++pos) {
    const unsigned prev = order[pos - 1];
    const unsigned a = order[pos];
    const unsigned b = order[pos + 1];
    const bool has_next = (pos + 2 < order.size());
    const unsigned next = has_next ? order[pos + 2] : 0;

    const int64_t before = uclaEdgeCost(cells[prev], cells[a])
                           + uclaEdgeCost(cells[a], cells[b])
                           + (has_next ? uclaEdgeCost(cells[b], cells[next]) : 0);
    const int64_t after = uclaEdgeCost(cells[prev], cells[b])
                          + uclaEdgeCost(cells[b], cells[a])
                          + (has_next ? uclaEdgeCost(cells[a], cells[next]) : 0);

    if (after < before
        || (after == before && cells[b].name < cells[a].name)) {
      std::swap(order[pos], order[pos + 1]);
      improved = true;
    }
  }
  return improved;
}

bool uclaWorstEdgeSegmentRelocateImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 4 || max_seg_len < 2) {
    return false;
  }

  const std::size_t j_limit = end_fixed ? (n - 1) : n;
  if (j_limit < 4) {
    return false;
  }

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost
        = uclaEdgeCost(cells[order[i]], cells[order[i + 1]]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  const std::size_t from = worst_i + 1;
  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (from == 0 || from > last_movable) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_to = 0;
  std::size_t best_after = 0;
  std::string_view best_name;

  const std::size_t seg_limit = std::min(max_seg_len, last_movable - from + 1);
  for (std::size_t seg_len = 2; seg_len <= seg_limit; ++seg_len) {
    if (uclaTimeExpired(deadline)) {
      break;
    }

    const std::size_t to = from + seg_len - 1;
    const unsigned prev = order[from - 1];
    const unsigned first = order[from];
    const unsigned last = order[to];
    const bool has_next = (to + 1 < n);
    const unsigned next = has_next ? order[to + 1] : 0;

    const int64_t remove_before = uclaEdgeCost(cells[prev], cells[first]);
    const int64_t remove_after
        = has_next ? uclaEdgeCost(cells[last], cells[next]) : 0;
    const int64_t remove_join
        = has_next ? uclaEdgeCost(cells[prev], cells[next]) : 0;
    const int64_t remove_delta = remove_join - remove_before - remove_after;

    for (std::size_t after = 0; after <= last_movable; ++after) {
      if ((after & 255U) == 0U && uclaTimeExpired(deadline)) {
        break;
      }

      // Disallow inserting inside the segment or immediately before it (no-op).
      if (after >= from - 1 && after <= to) {
        continue;
      }

      const unsigned ins_prev = order[after];
      const bool has_ins_next = (after + 1 < n);
      const unsigned ins_next = has_ins_next ? order[after + 1] : 0;

      const int64_t insert_before
          = has_ins_next ? uclaEdgeCost(cells[ins_prev], cells[ins_next]) : 0;
      const int64_t insert_after
          = uclaEdgeCost(cells[ins_prev], cells[first])
            + (has_ins_next ? uclaEdgeCost(cells[last], cells[ins_next]) : 0);
      const int64_t insert_delta = insert_after - insert_before;

      const int64_t delta = remove_delta + insert_delta;
      const std::string_view key = cells[first].name;
      if (delta < best_delta
          || (delta == best_delta && found
              && (key < best_name
                  || (key == best_name
                      && (to < best_to
                          || (to == best_to && after < best_after)))))) {
        found = true;
        best_delta = delta;
        best_to = to;
        best_after = after;
        best_name = key;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::vector<unsigned> segment(
      order.begin() + static_cast<std::ptrdiff_t>(from),
      order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));
  order.erase(order.begin() + static_cast<std::ptrdiff_t>(from),
              order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));

  std::size_t after = best_after;
  if (after > best_to) {
    after -= segment.size();
  }
  const std::size_t insert_index = after + 1;
  order.insert(order.begin() + static_cast<std::ptrdiff_t>(insert_index),
               segment.begin(),
               segment.end());
  return true;
}

bool uclaWorstEdgeSegmentSwapImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 4 || max_seg_len == 0) {
    return false;
  }

  // Relocation needs a successor node, so keep the "last movable" position at
  // n-2 regardless of whether endpoints are fixed.
  const std::size_t last_movable = n - 2;
  if (last_movable < 2) {
    return false;
  }

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 1; i + 1 <= last_movable; ++i) {
    const int64_t cost
        = uclaEdgeCost(cells[order[i]], cells[order[i + 1]]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  const std::size_t b = worst_i + 1;
  if (b < 2 || b > last_movable) {
    return false;
  }

  const std::size_t max_left = std::min(max_seg_len, b - 1);
  const std::size_t max_right = std::min(max_seg_len, last_movable - b + 1);
  if (max_left == 0 || max_right == 0) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_len1 = 0;
  std::size_t best_len2 = 0;
  std::string_view best_key;

  const unsigned x = order[b - 1];
  const unsigned bnode = order[b];

  for (std::size_t len1 = 1; len1 <= max_left; ++len1) {
    if (uclaTimeExpired(deadline)) {
      break;
    }
    const std::size_t a = b - len1;
    if (a < 1) {
      continue;
    }
    const unsigned p = order[a - 1];
    const unsigned anode = order[a];

    for (std::size_t len2 = 1; len2 <= max_right; ++len2) {
      if ((len2 & 1023U) == 0U && uclaTimeExpired(deadline)) {
        break;
      }

      const std::size_t c = b + len2 - 1;
      if (c < b || c > last_movable) {
        continue;
      }

      const unsigned cnode = order[c];
      const bool has_n = (c + 1 < n);
      const unsigned nnode = has_n ? order[c + 1] : 0;

      const int64_t before = uclaEdgeCost(cells[p], cells[anode])
                             + uclaEdgeCost(cells[x], cells[bnode])
                             + (has_n ? uclaEdgeCost(cells[cnode], cells[nnode])
                                      : 0);
      const int64_t after = uclaEdgeCost(cells[p], cells[bnode])
                            + uclaEdgeCost(cells[cnode], cells[anode])
                            + (has_n ? uclaEdgeCost(cells[x], cells[nnode]) : 0);
      const int64_t delta = after - before;

      const std::string_view key = cells[bnode].name;
      if (delta < best_delta
          || (delta == best_delta && found
              && (key < best_key
                  || (key == best_key
                      && (len1 < best_len1
                          || (len1 == best_len1 && len2 < best_len2)))))) {
        found = true;
        best_delta = delta;
        best_len1 = len1;
        best_len2 = len2;
        best_key = key;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  const std::size_t a = b - best_len1;
  const std::size_t c = b + best_len2 - 1;
  std::rotate(order.begin() + static_cast<std::ptrdiff_t>(a),
              order.begin() + static_cast<std::ptrdiff_t>(b),
              order.begin() + static_cast<std::ptrdiff_t>(c + 1));
  return true;
}

bool uclaWorstEdgeTwoOptImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  // j indexes the start of the second edge (c->d), so j+1 must be valid.
  const std::size_t j_limit = n - 1;

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost
        = uclaEdgeCost(cells[order[i]], cells[order[i + 1]]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  // Prefix sums for forward edges and reversed-adjacent edges.
  std::vector<int64_t> forward_prefix(n, 0);
  std::vector<int64_t> rev_prefix(n, 0);
  for (std::size_t i = 1; i < n; ++i) {
    forward_prefix[i] = forward_prefix[i - 1]
                        + uclaEdgeCost(cells[order[i - 1]], cells[order[i]]);
    rev_prefix[i] = rev_prefix[i - 1]
                    + uclaEdgeCost(cells[order[i]], cells[order[i - 1]]);
  }

  const auto segment_forward = [&](std::size_t from,
                                   std::size_t to) -> int64_t {
    return forward_prefix[to] - forward_prefix[from];
  };
  const auto segment_reversed = [&](std::size_t from,
                                    std::size_t to) -> int64_t {
    return rev_prefix[to] - rev_prefix[from];
  };

  const std::size_t i = worst_i;
  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_j = 0;

  const std::size_t max_j = end_fixed ? j_limit : j_limit;
  for (std::size_t j = i + 2; j < max_j; ++j) {
    if ((j & 255U) == 0U && uclaTimeExpired(deadline)) {
      break;
    }

    const unsigned a = order[i];
    const unsigned b = order[i + 1];
    const unsigned c = order[j];
    const unsigned d = order[j + 1];

    const int64_t old_inside = segment_forward(i + 1, j);
    const int64_t new_inside = segment_reversed(i + 1, j);

    const int64_t old_edges
        = uclaEdgeCost(cells[a], cells[b]) + uclaEdgeCost(cells[c], cells[d]);
    const int64_t new_edges
        = uclaEdgeCost(cells[a], cells[c]) + uclaEdgeCost(cells[b], cells[d]);

    const int64_t before = old_edges + old_inside;
    const int64_t after = new_edges + new_inside;
    const int64_t delta = after - before;

    const std::string_view key_c = cells[c].name;
    const std::string_view best_key = found ? cells[order[best_j]].name : "";
    if (delta < best_delta
        || (delta == best_delta && found
            && (key_c < best_key
                || (key_c == best_key
                    && cells[a].name < cells[order[i]].name)))) {
      found = true;
      best_delta = delta;
      best_j = j;
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::reverse(order.begin() + static_cast<std::ptrdiff_t>(i + 1),
               order.begin() + static_cast<std::ptrdiff_t>(best_j + 1));
  return true;
}

bool uclaTwoOptBestImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  const std::size_t last_pos = end_fixed ? (n - 2) : (n - 1);
  if (last_pos < 2) {
    return false;
  }

  // Prefix sums for forward edges and reversed-adjacent edges.
  std::vector<int64_t> forward_prefix(n, 0);
  std::vector<int64_t> rev_prefix(n, 0);
  for (std::size_t i = 1; i < n; ++i) {
    forward_prefix[i] = forward_prefix[i - 1]
                        + uclaEdgeCost(cells[order[i - 1]], cells[order[i]]);
    rev_prefix[i] = rev_prefix[i - 1]
                    + uclaEdgeCost(cells[order[i]], cells[order[i - 1]]);
  }

  const auto segment_forward = [&](std::size_t from,
                                   std::size_t to) -> int64_t {
    return forward_prefix[to] - forward_prefix[from];
  };
  const auto segment_reversed = [&](std::size_t from,
                                    std::size_t to) -> int64_t {
    return rev_prefix[to] - rev_prefix[from];
  };

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_i = 0;
  std::size_t best_j = 0;

  for (std::size_t i = 0; i + 2 <= last_pos; ++i) {
    if ((i & 63U) == 0U && uclaTimeExpired(deadline)) {
      break;
    }
    const unsigned a = order[i];
    const unsigned b = order[i + 1];

    for (std::size_t j = i + 2; j <= last_pos; ++j) {
      if ((j & 255U) == 0U && uclaTimeExpired(deadline)) {
        break;
      }

      const unsigned c = order[j];
      const unsigned d = order[j + 1];

      const int64_t old_inside = segment_forward(i + 1, j);
      const int64_t new_inside = segment_reversed(i + 1, j);

      const int64_t old_edges = uclaEdgeCost(cells[a], cells[b])
                                + uclaEdgeCost(cells[c], cells[d]);
      const int64_t new_edges = uclaEdgeCost(cells[a], cells[c])
                                + uclaEdgeCost(cells[b], cells[d]);

      const int64_t before = old_edges + old_inside;
      const int64_t after = new_edges + new_inside;
      const int64_t delta = after - before;

      if (delta < best_delta) {
        found = true;
        best_delta = delta;
        best_i = i;
        best_j = j;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::reverse(order.begin() + static_cast<std::ptrdiff_t>(best_i + 1),
               order.begin() + static_cast<std::ptrdiff_t>(best_j + 1));
  return true;
}

bool uclaOrOptBestImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 4 || max_seg_len == 0) {
    return false;
  }

  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (last_movable < 2) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_from = 0;
  std::size_t best_to = 0;
  std::size_t best_after = 0;
  std::string_view best_key;

  const std::size_t seg_limit = std::min(max_seg_len, last_movable);
  for (std::size_t seg_len = 1; seg_len <= seg_limit; ++seg_len) {
    if (uclaTimeExpired(deadline)) {
      break;
    }

    for (std::size_t from = 1; from + seg_len - 1 <= last_movable; ++from) {
      const std::size_t to = from + seg_len - 1;
      const unsigned prev = order[from - 1];
      const unsigned first = order[from];
      const unsigned last = order[to];
      const unsigned next = order[to + 1];

      const int64_t remove_before = uclaEdgeCost(cells[prev], cells[first]);
      const int64_t remove_after = uclaEdgeCost(cells[last], cells[next]);
      const int64_t remove_join = uclaEdgeCost(cells[prev], cells[next]);
      const int64_t remove_delta = remove_join - remove_before - remove_after;

      for (std::size_t after = 0; after <= last_movable; ++after) {
        if ((after & 255U) == 0U && uclaTimeExpired(deadline)) {
          break;
        }

        // Disallow inserting inside the segment or immediately before it (no-op).
        if (after >= from - 1 && after <= to) {
          continue;
        }

        const unsigned ins_prev = order[after];
        const unsigned ins_next = order[after + 1];

        const int64_t insert_before = uclaEdgeCost(cells[ins_prev], cells[ins_next]);
        const int64_t insert_after = uclaEdgeCost(cells[ins_prev], cells[first])
                                     + uclaEdgeCost(cells[last], cells[ins_next]);
        const int64_t insert_delta = insert_after - insert_before;

        const int64_t delta = remove_delta + insert_delta;
        const std::string_view key = cells[first].name;
        bool take = false;
        if (delta < best_delta) {
          take = true;
        } else if (delta == best_delta && found) {
          if (key < best_key) {
            take = true;
          } else if (key == best_key) {
            const std::size_t best_seg_len = best_to - best_from + 1;
            if (seg_len < best_seg_len) {
              take = true;
            } else if (seg_len == best_seg_len) {
              if (from < best_from || (from == best_from && after < best_after)) {
                take = true;
              }
            }
          }
        }

        if (take) {
          found = true;
          best_delta = delta;
          best_from = from;
          best_to = to;
          best_after = after;
          best_key = key;
        }
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::vector<unsigned> segment(
      order.begin() + static_cast<std::ptrdiff_t>(best_from),
      order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));
  order.erase(order.begin() + static_cast<std::ptrdiff_t>(best_from),
              order.begin() + static_cast<std::ptrdiff_t>(best_to + 1));

  std::size_t after = best_after;
  if (after > best_to) {
    after -= segment.size();
  }
  const std::size_t insert_index = after + 1;
  order.insert(order.begin() + static_cast<std::ptrdiff_t>(insert_index),
               segment.begin(),
               segment.end());
  return true;
}

bool uclaBestSwapImprove(const std::vector<FlatScanChain::CellData>& cells,
                         std::vector<unsigned>& order,
                         bool end_fixed)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  bool found = false;
  int64_t best_delta = 0;
  std::size_t best_i = 0;
  std::size_t best_j = 0;
  const std::size_t j_limit = end_fixed ? (n - 1) : n;

  for (std::size_t i = 1; i + 1 < n; ++i) {
    for (std::size_t j = i + 1; j < j_limit; ++j) {
      const unsigned ai = order[i];
      const unsigned aj = order[j];

      const bool has_pi = (i > 0);
      const bool has_ni = (i + 1 < n);
      const bool has_pj = (j > 0);
      const bool has_nj = (j + 1 < n);

      const unsigned pi = has_pi ? order[i - 1] : 0;
      const unsigned ni = has_ni ? order[i + 1] : 0;
      const unsigned pj = has_pj ? order[j - 1] : 0;
      const unsigned nj = has_nj ? order[j + 1] : 0;

      int64_t before = 0;
      int64_t after = 0;

      if (has_pi) {
        before += uclaEdgeCost(cells[pi], cells[ai]);
        after += uclaEdgeCost(cells[pi], cells[aj]);
      }
      if (has_ni) {
        if (i + 1 == j) {
          before += uclaEdgeCost(cells[ai], cells[aj]);
          after += uclaEdgeCost(cells[aj], cells[ai]);
        } else {
          before += uclaEdgeCost(cells[ai], cells[ni]);
          after += uclaEdgeCost(cells[aj], cells[ni]);
        }
      }

      if (has_pj && j - 1 != i) {
        before += uclaEdgeCost(cells[pj], cells[aj]);
        after += uclaEdgeCost(cells[pj], cells[ai]);
      }
      if (has_nj) {
        before += uclaEdgeCost(cells[aj], cells[nj]);
        after += uclaEdgeCost(cells[ai], cells[nj]);
      }

      const int64_t delta = after - before;
      const bool take
          = delta < best_delta
            || (delta == best_delta && found
                && (cells[aj].name < cells[order[best_j]].name
                    || (cells[aj].name == cells[order[best_j]].name
                        && cells[ai].name < cells[order[best_i]].name)));
      if (take) {
        found = true;
        best_delta = delta;
        best_i = i;
        best_j = j;
      }
    }
  }

  if (!found || best_delta >= 0) {
    return false;
  }

  std::swap(order[best_i], order[best_j]);
  return true;
}

void uclaRefinePathL1(const std::vector<FlatScanChain::CellData>& cells,
                      std::vector<unsigned>& order,
                      bool end_fixed,
                      const std::optional<std::chrono::steady_clock::time_point>&
                          deadline)
{
  // Short deterministic local descent; bounded by the overall deadline.
  for (int pass = 0; pass < 200; ++pass) {
    if (uclaTimeExpired(deadline)) {
      break;
    }
    bool improved = uclaAdjacentSwapImprove(cells, order, end_fixed);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaWorstEdgeTwoOptImprove(cells, order, end_fixed, deadline);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaWorstEdgeSegmentRelocateImprove(
        cells, order, end_fixed, deadline, 256);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaWorstEdgeSegmentSwapImprove(cells, order, end_fixed, deadline, 256);
    if (!improved) {
      break;
    }
  }
}

void uclaRefinePathL1Intensive(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  // More thorough local descent; still bounded by the overall deadline.
  for (int pass = 0; pass < 2000; ++pass) {
    if (uclaTimeExpired(deadline)) {
      break;
    }

    bool improved = false;
    improved |= uclaTwoOptBestImprove(cells, order, end_fixed, deadline);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaOrOptBestImprove(cells, order, end_fixed, deadline, 3);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaAdjacentSwapImprove(cells, order, end_fixed);
    if (uclaTimeExpired(deadline)) {
      break;
    }

    if (!improved) {
      break;
    }
  }
}

bool uclaWorstEdgeTwoOptScoreImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  const std::size_t n = order.size();
  if (n < 4) {
    return false;
  }

  const std::size_t j_limit = n - 1;
  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost = uclaEdgeCost(cells[order[i]], cells[order[i + 1]]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  std::vector<int64_t> forward_prefix(n, 0);
  std::vector<int64_t> rev_prefix(n, 0);
  for (std::size_t i = 1; i < n; ++i) {
    forward_prefix[i] = forward_prefix[i - 1]
                        + uclaEdgeCost(cells[order[i - 1]], cells[order[i]]);
    rev_prefix[i] = rev_prefix[i - 1]
                    + uclaEdgeCost(cells[order[i]], cells[order[i - 1]]);
  }

  const auto segment_forward = [&](std::size_t from,
                                   std::size_t to) -> int64_t {
    return forward_prefix[to] - forward_prefix[from];
  };
  const auto segment_reversed = [&](std::size_t from,
                                    std::size_t to) -> int64_t {
    return rev_prefix[to] - rev_prefix[from];
  };

  const UclaPathScore base_score = uclaEvaluatePath(cells, order);
  UclaPathScore best_score = base_score;
  std::vector<unsigned> best_path;
  bool found = false;

  const std::size_t i = worst_i;
  const std::size_t max_j = end_fixed ? j_limit : j_limit;
  for (std::size_t j = i + 2; j < max_j; ++j) {
    if ((j & 255U) == 0U && uclaTimeExpired(deadline)) {
      break;
    }

    const unsigned a = order[i];
    const unsigned b = order[i + 1];
    const unsigned c = order[j];
    const unsigned d = order[j + 1];

    const int64_t old_inside = segment_forward(i + 1, j);
    const int64_t new_inside = segment_reversed(i + 1, j);
    const int64_t old_edges
        = uclaEdgeCost(cells[a], cells[b]) + uclaEdgeCost(cells[c], cells[d]);
    const int64_t new_edges
        = uclaEdgeCost(cells[a], cells[c]) + uclaEdgeCost(cells[b], cells[d]);
    const int64_t delta = (new_edges + new_inside) - (old_edges + old_inside);
    if (delta > 0) {
      continue;
    }

    std::vector<unsigned> cand = order;
    std::reverse(cand.begin() + static_cast<std::ptrdiff_t>(i + 1),
                 cand.begin() + static_cast<std::ptrdiff_t>(j + 1));
    const UclaPathScore cand_score = uclaEvaluatePath(cells, cand);
    const bool better_metric = uclaMetricLess(cand_score, best_score);
    const bool same_metric = cand_score.cost == best_score.cost
                             && cand_score.max_edge == best_score.max_edge
                             && cand_score.tail_sum == best_score.tail_sum;
    if (better_metric || (same_metric && found && std::lexicographical_compare(cand.begin(), cand.end(), best_path.begin(), best_path.end()))) {
      found = true;
      best_score = cand_score;
      best_path = std::move(cand);
    }
  }

  if (!found) {
    return false;
  }

  order.swap(best_path);
  return true;
}

bool uclaWorstEdgeSegmentRelocateScoreImprove(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::size_t max_seg_len)
{
  const std::size_t n = order.size();
  if (n < 3 || max_seg_len == 0) {
    return false;
  }

  const std::size_t j_limit = end_fixed ? (n - 1) : n;
  if (j_limit < 3) {
    return false;
  }

  std::size_t worst_i = 0;
  int64_t worst_cost = -1;
  bool have_worst = false;
  for (std::size_t i = 0; i + 2 < j_limit; ++i) {
    const int64_t cost = uclaEdgeCost(cells[order[i]], cells[order[i + 1]]);
    if (!have_worst || cost > worst_cost) {
      have_worst = true;
      worst_cost = cost;
      worst_i = i;
    }
  }
  if (!have_worst) {
    return false;
  }

  const std::size_t from = worst_i + 1;
  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (from == 0 || from > last_movable) {
    return false;
  }

  const UclaPathScore base_score = uclaEvaluatePath(cells, order);
  UclaPathScore best_score = base_score;
  std::vector<unsigned> best_path;
  bool found = false;

  const std::size_t seg_limit = std::min(max_seg_len, last_movable - from + 1);
  for (std::size_t seg_len = 1; seg_len <= seg_limit; ++seg_len) {
    if (uclaTimeExpired(deadline)) {
      break;
    }

    const std::size_t to = from + seg_len - 1;
    const unsigned prev = order[from - 1];
    const unsigned first = order[from];
    const unsigned last = order[to];
    const bool has_next = (to + 1 < n);
    const unsigned next = has_next ? order[to + 1] : 0;

    const int64_t remove_before = uclaEdgeCost(cells[prev], cells[first]);
    const int64_t remove_after
        = has_next ? uclaEdgeCost(cells[last], cells[next]) : 0;
    const int64_t remove_join
        = has_next ? uclaEdgeCost(cells[prev], cells[next]) : 0;
    const int64_t remove_delta = remove_join - remove_before - remove_after;

    for (std::size_t after = 0; after <= last_movable; ++after) {
      if ((after & 255U) == 0U && uclaTimeExpired(deadline)) {
        break;
      }
      if (after >= from - 1 && after <= to) {
        continue;
      }

      const unsigned ins_prev = order[after];
      const bool has_ins_next = (after + 1 < n);
      const unsigned ins_next = has_ins_next ? order[after + 1] : 0;

      const int64_t insert_before
          = has_ins_next ? uclaEdgeCost(cells[ins_prev], cells[ins_next]) : 0;
      const int64_t insert_after
          = uclaEdgeCost(cells[ins_prev], cells[first])
            + (has_ins_next ? uclaEdgeCost(cells[last], cells[ins_next]) : 0);
      const int64_t delta = remove_delta + (insert_after - insert_before);
      if (delta > 0) {
        continue;
      }

      std::vector<unsigned> cand = order;
      std::vector<unsigned> segment(
          cand.begin() + static_cast<std::ptrdiff_t>(from),
          cand.begin() + static_cast<std::ptrdiff_t>(to + 1));
      cand.erase(cand.begin() + static_cast<std::ptrdiff_t>(from),
                 cand.begin() + static_cast<std::ptrdiff_t>(to + 1));

      std::size_t after_adj = after;
      if (after_adj > to) {
        after_adj -= segment.size();
      }
      const std::size_t insert_index = after_adj + 1;
      cand.insert(cand.begin() + static_cast<std::ptrdiff_t>(insert_index),
                  segment.begin(),
                  segment.end());

      const UclaPathScore cand_score = uclaEvaluatePath(cells, cand);
      const bool better_metric = uclaMetricLess(cand_score, best_score);
      const bool same_metric = cand_score.cost == best_score.cost
                               && cand_score.max_edge == best_score.max_edge
                               && cand_score.tail_sum == best_score.tail_sum;
      if (better_metric || (same_metric && found && std::lexicographical_compare(cand.begin(), cand.end(), best_path.begin(), best_path.end()))) {
        found = true;
        best_score = cand_score;
        best_path = std::move(cand);
      }
    }
  }

  if (!found) {
    return false;
  }

  order.swap(best_path);
  return true;
}

void uclaRefinePathJumpPlateau(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  for (int pass = 0; pass < 6; ++pass) {
    if (uclaTimeExpired(deadline)) {
      break;
    }
    bool improved
        = uclaWorstEdgeTwoOptScoreImprove(cells, order, end_fixed, deadline);
    if (uclaTimeExpired(deadline)) {
      break;
    }
    improved |= uclaWorstEdgeSegmentRelocateScoreImprove(
        cells, order, end_fixed, deadline, 3);
    if (!improved) {
      break;
    }
  }
}

int64_t uclaTailSumFromEdges(const std::array<int64_t, 4>& edges)
{
  int64_t sum = 0;
  for (const int64_t edge : edges) {
    sum += edge;
  }
  return sum;
}

struct UclaExactTailState
{
  int64_t cost = std::numeric_limits<int64_t>::max();
  std::array<int64_t, 4> worst_edges{};
  bool valid = false;
};

UclaPathScore uclaExactTailStateToScore(const UclaExactTailState& state)
{
  UclaPathScore score;
  score.cost = state.cost;
  score.max_edge = state.worst_edges[0];
  score.tail_sum = uclaTailSumFromEdges(state.worst_edges);
  return score;
}

bool uclaExactTailStateLess(const UclaExactTailState& a,
                            const UclaExactTailState& b)
{
  if (a.valid != b.valid) {
    return a.valid && !b.valid;
  }
  if (!a.valid) {
    return false;
  }

  const UclaPathScore score_a = uclaExactTailStateToScore(a);
  const UclaPathScore score_b = uclaExactTailStateToScore(b);
  if (uclaMetricLess(score_a, score_b)) {
    return true;
  }
  if (uclaMetricLess(score_b, score_a)) {
    return false;
  }

  for (std::size_t i = 0; i < a.worst_edges.size(); ++i) {
    if (a.worst_edges[i] != b.worst_edges[i]) {
      return a.worst_edges[i] < b.worst_edges[i];
    }
  }
  return false;
}

UclaExactTailState uclaAppendTailEdge(UclaExactTailState state, int64_t edge)
{
  if (!state.valid) {
    return state;
  }

  state.cost += edge;
  for (std::size_t pos = 0; pos < state.worst_edges.size(); ++pos) {
    if (edge <= state.worst_edges[pos]) {
      continue;
    }
    for (std::size_t shift = state.worst_edges.size(); shift > pos + 1; --shift) {
      state.worst_edges[shift - 1] = state.worst_edges[shift - 2];
    }
    state.worst_edges[pos] = edge;
    break;
  }
  return state;
}

std::vector<std::size_t> uclaCollectWorstEdgePositions(
    const std::vector<FlatScanChain::CellData>& cells,
    const std::vector<unsigned>& order,
    std::size_t max_count)
{
  struct WorstEdgeEntry
  {
    int64_t edge = -1;
    std::size_t pos = 0;
  };

  std::vector<WorstEdgeEntry> worst;
  worst.reserve(max_count);

  for (std::size_t pos = 0; pos + 1 < order.size(); ++pos) {
    const int64_t edge = uclaEdgeCost(cells[order[pos]], cells[order[pos + 1]]);
    WorstEdgeEntry entry{edge, pos};

    auto it = worst.begin();
    while (it != worst.end()
           && (it->edge > edge || (it->edge == edge && it->pos < pos))) {
      ++it;
    }

    if (it != worst.end() || worst.size() < max_count) {
      worst.insert(it, entry);
      if (worst.size() > max_count) {
        worst.pop_back();
      }
    }
  }

  std::vector<std::size_t> out;
  out.reserve(worst.size());
  for (const auto& entry : worst) {
    out.push_back(entry.pos);
  }
  return out;
}

std::size_t uclaClampWindowStart(std::size_t center,
                                 std::size_t window_size,
                                 std::size_t last_movable)
{
  if (window_size >= last_movable) {
    return 1;
  }

  const std::size_t max_start = last_movable + 1 - window_size;
  std::size_t start = 1;
  if (center > window_size / 2) {
    start = center - window_size / 2;
  }
  if (start > max_start) {
    start = max_start;
  }
  return std::max<std::size_t>(1, start);
}

bool uclaExactWindowImprove(const std::vector<FlatScanChain::CellData>& cells,
                            std::vector<unsigned>& order,
                            UclaPathScore& current_score,
                            std::size_t from,
                            std::size_t to,
                            const std::optional<std::chrono::steady_clock::time_point>&
                                deadline)
{
  constexpr std::size_t kMaxExactWindow = 14;

  if (from == 0 || from > to || to + 1 >= order.size()) {
    return false;
  }

  const std::size_t count = to - from + 1;
  if (count < 2 || count > kMaxExactWindow) {
    return false;
  }
  if (uclaTimeExpired(deadline)) {
    return false;
  }

  const unsigned left_anchor = order[from - 1];
  const unsigned right_anchor = order[to + 1];

  std::vector<unsigned> window_nodes(
      order.begin() + static_cast<std::ptrdiff_t>(from),
      order.begin() + static_cast<std::ptrdiff_t>(to + 1));

  std::vector<unsigned> current_local_path;
  current_local_path.reserve(count + 2);
  current_local_path.push_back(left_anchor);
  current_local_path.insert(
      current_local_path.end(), window_nodes.begin(), window_nodes.end());
  current_local_path.push_back(right_anchor);
  const UclaPathScore current_local_score
      = uclaEvaluatePath(cells, current_local_path);

  const std::size_t states = static_cast<std::size_t>(1U) << count;
  std::vector<int64_t> start_cost(count, 0);
  std::vector<int64_t> end_cost(count, 0);
  std::vector<int64_t> edge_costs(count * count, 0);
  for (std::size_t i = 0; i < count; ++i) {
    start_cost[i] = uclaEdgeCost(cells[left_anchor], cells[window_nodes[i]]);
    end_cost[i] = uclaEdgeCost(cells[window_nodes[i]], cells[right_anchor]);
    for (std::size_t j = 0; j < count; ++j) {
      edge_costs[i * count + j]
          = uclaEdgeCost(cells[window_nodes[i]], cells[window_nodes[j]]);
    }
  }

  std::vector<UclaExactTailState> dp(states * count);
  std::vector<int16_t> prev_last(states * count, -1);
  std::vector<uint16_t> prev_mask(states * count, 0);

  const auto index_of = [&](std::size_t mask, std::size_t last) {
    return mask * count + last;
  };

  UclaExactTailState empty;
  empty.valid = true;
  empty.cost = 0;
  empty.worst_edges.fill(0);

  for (std::size_t last = 0; last < count; ++last) {
    const std::size_t mask = static_cast<std::size_t>(1U) << last;
    dp[index_of(mask, last)] = uclaAppendTailEdge(empty, start_cost[last]);
  }

  for (std::size_t mask = 1; mask < states; ++mask) {
    if ((mask & 127U) == 0U && uclaTimeExpired(deadline)) {
      break;
    }

    for (std::size_t last = 0; last < count; ++last) {
      if ((mask & (static_cast<std::size_t>(1U) << last)) == 0U) {
        continue;
      }
      const UclaExactTailState cur = dp[index_of(mask, last)];
      if (!cur.valid) {
        continue;
      }

      for (std::size_t next = 0; next < count; ++next) {
        const std::size_t next_bit = static_cast<std::size_t>(1U) << next;
        if ((mask & next_bit) != 0U) {
          continue;
        }

        const std::size_t next_mask = mask | next_bit;
        const UclaExactTailState cand
            = uclaAppendTailEdge(cur, edge_costs[last * count + next]);
        const std::size_t next_index = index_of(next_mask, next);
        if (uclaExactTailStateLess(cand, dp[next_index])) {
          dp[next_index] = cand;
          prev_last[next_index] = static_cast<int16_t>(last);
          prev_mask[next_index] = static_cast<uint16_t>(mask);
        }
      }
    }
  }

  const std::size_t full_mask = states - 1;
  UclaExactTailState best_state;
  int best_last = -1;
  for (std::size_t last = 0; last < count; ++last) {
    const UclaExactTailState cur = dp[index_of(full_mask, last)];
    if (!cur.valid) {
      continue;
    }

    const UclaExactTailState cand = uclaAppendTailEdge(cur, end_cost[last]);
    if (uclaExactTailStateLess(cand, best_state)) {
      best_state = cand;
      best_last = static_cast<int>(last);
    }
  }

  if (!best_state.valid || best_last < 0) {
    return false;
  }

  std::vector<unsigned> best_segment(count, 0);
  std::size_t mask = full_mask;
  int last = best_last;
  for (std::size_t pos = count; pos > 0; --pos) {
    if (last < 0) {
      return false;
    }
    best_segment[pos - 1] = window_nodes[static_cast<std::size_t>(last)];
    const std::size_t idx = index_of(mask, static_cast<std::size_t>(last));
    const std::size_t next_mask = prev_mask[idx];
    last = prev_last[idx];
    mask = next_mask;
  }

  std::vector<unsigned> best_local_path;
  best_local_path.reserve(count + 2);
  best_local_path.push_back(left_anchor);
  best_local_path.insert(
      best_local_path.end(), best_segment.begin(), best_segment.end());
  best_local_path.push_back(right_anchor);
  const UclaPathScore best_local_score = uclaEvaluatePath(cells, best_local_path);
  if (!uclaScoreLess(best_local_score,
                     best_local_path,
                     current_local_score,
                     current_local_path)) {
    return false;
  }

  std::vector<unsigned> candidate = order;
  std::copy(best_segment.begin(),
            best_segment.end(),
            candidate.begin() + static_cast<std::ptrdiff_t>(from));
  const UclaPathScore candidate_score = uclaEvaluatePath(cells, candidate);
  const bool global_better
      = uclaScoreLess(candidate_score, candidate, current_score, order);
  const bool global_same_metric = candidate_score.cost == current_score.cost
                                  && candidate_score.max_edge == current_score.max_edge
                                  && candidate_score.tail_sum == current_score.tail_sum;
  if (!global_better && !global_same_metric) {
    return false;
  }

  order.swap(candidate);
  current_score = candidate_score;
  return true;
}

bool uclaRefinePathExactWindows(
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& order,
    bool end_fixed,
    std::mt19937_64& rng,
    const std::optional<std::chrono::steady_clock::time_point>& deadline)
{
  const std::size_t n = order.size();
  if (n < 5) {
    return false;
  }

  const std::size_t last_movable = end_fixed ? (n - 2) : (n - 1);
  if (last_movable < 2) {
    return false;
  }

  std::vector<std::size_t> window_sizes;
  window_sizes.reserve(3);
  auto push_window_size = [&](std::size_t size) {
    if (size < 2 || size > last_movable) {
      return;
    }
    if (std::find(window_sizes.begin(), window_sizes.end(), size)
        == window_sizes.end()) {
      window_sizes.push_back(size);
    }
  };
  push_window_size(std::min<std::size_t>(14, last_movable));
  push_window_size(std::min<std::size_t>(12, last_movable));
  push_window_size(std::min<std::size_t>(10, last_movable));
  if (window_sizes.empty()) {
    return false;
  }

  UclaPathScore current_score = uclaEvaluatePath(cells, order);
  bool any_improved = false;

  constexpr std::size_t kWorstEdgeCount = 8;
  constexpr std::size_t kRandomWindowsPerSize = 1;
  for (int pass = 0; pass < 4; ++pass) {
    if (uclaTimeExpired(deadline)) {
      break;
    }

    bool pass_improved = false;
    const std::vector<std::size_t> worst_positions
        = uclaCollectWorstEdgePositions(cells, order, kWorstEdgeCount);
    std::vector<std::pair<std::size_t, std::size_t>> attempted_windows;
    attempted_windows.reserve(window_sizes.size() * (worst_positions.size() + 2));

    auto try_window = [&](std::size_t start, std::size_t size) {
      if (size < 2 || start == 0 || start + size - 1 > last_movable) {
        return;
      }
      const std::pair<std::size_t, std::size_t> key{start, size};
      if (std::find(attempted_windows.begin(), attempted_windows.end(), key)
          != attempted_windows.end()) {
        return;
      }
      attempted_windows.push_back(key);
      if (uclaExactWindowImprove(
              cells, order, current_score, start, start + size - 1, deadline)) {
        pass_improved = true;
      }
    };

    for (const std::size_t size : window_sizes) {
      for (const std::size_t edge_pos : worst_positions) {
        if (uclaTimeExpired(deadline)) {
          break;
        }
        const std::size_t center = std::min(last_movable, edge_pos + 1);
        const std::size_t start
            = uclaClampWindowStart(center, size, last_movable);
        try_window(start, size);
      }
      if (uclaTimeExpired(deadline)) {
        break;
      }

      const std::size_t span = last_movable + 1 - size;
      if (span == 0) {
        try_window(1, size);
      } else {
        for (std::size_t sample = 0; sample < kRandomWindowsPerSize; ++sample) {
          const std::size_t start
              = 1 + static_cast<std::size_t>(rng() % span);
          try_window(start, size);
        }
      }
    }

    if (!pass_improved) {
      break;
    }
    any_improved = true;

    if (!uclaTimeExpired(deadline)) {
      uclaRefinePathJumpPlateau(cells, order, end_fixed, deadline);
      current_score = uclaEvaluatePath(cells, order);
    }
  }

  return any_improved;
}

struct MultiStartResult
{
  int64_t best_cost = std::numeric_limits<int64_t>::max();
  std::vector<unsigned> best_path;
};

struct EliteEntry
{
  int64_t cost = std::numeric_limits<int64_t>::max();
  std::vector<unsigned> path;
};

bool lexLess(const std::vector<unsigned>& a, const std::vector<unsigned>& b)
{
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

void updateElitePool(std::vector<EliteEntry>& elite,
                     std::size_t max_elite,
                     int64_t cost,
                     const std::vector<unsigned>& path)
{
  if (max_elite == 0) {
    return;
  }

  if (elite.size() >= max_elite) {
    const EliteEntry& worst = elite.back();
    if (cost > worst.cost || (cost == worst.cost && !lexLess(path, worst.path))) {
      return;
    }
  }

  EliteEntry e;
  e.cost = cost;
  e.path = path;
  elite.push_back(std::move(e));
  std::sort(elite.begin(),
            elite.end(),
            [](const EliteEntry& a, const EliteEntry& b) {
              if (a.cost != b.cost) {
                return a.cost < b.cost;
              }
              return lexLess(a.path, b.path);
            });
  if (elite.size() > max_elite) {
    elite.resize(max_elite);
  }
}

template <typename KeyFn>
std::vector<unsigned> preSortInterior(const std::vector<unsigned>& interior,
                                      KeyFn key_fn)
{
  std::vector<unsigned> out = interior;
  std::sort(out.begin(), out.end(), [&](unsigned a, unsigned b) {
    const auto ka = key_fn(a);
    const auto kb = key_fn(b);
    if (ka.first != kb.first) {
      return ka.first < kb.first;
    }
    if (ka.second != kb.second) {
      return ka.second < kb.second;
    }
    // Deterministic tiebreak.
    return a < b;
  });
  return out;
}

void buildPathFromInterior(unsigned begin_idx,
                           unsigned end_idx,
                           const std::vector<unsigned>& interior,
                           std::vector<unsigned>& out_path)
{
  out_path.clear();
  out_path.reserve(interior.size() + 2);
  out_path.push_back(begin_idx);
  out_path.insert(out_path.end(), interior.begin(), interior.end());
  out_path.push_back(end_idx);
}

void buildChunkShuffledInterior(const std::vector<unsigned>& base_interior,
                                std::size_t chunk_size,
                                std::mt19937_64& rng,
                                std::vector<unsigned>& out_interior)
{
  out_interior.clear();
  out_interior.reserve(base_interior.size());
  if (base_interior.empty()) {
    return;
  }
  if (chunk_size == 0) {
    out_interior = base_interior;
    return;
  }

  const std::size_t n = base_interior.size();
  const std::size_t n_chunks = (n + chunk_size - 1) / chunk_size;
  std::vector<std::size_t> chunk_ids;
  chunk_ids.reserve(n_chunks);
  for (std::size_t i = 0; i < n_chunks; ++i) {
    chunk_ids.push_back(i);
  }
  std::shuffle(chunk_ids.begin(), chunk_ids.end(), rng);

  for (const std::size_t cid : chunk_ids) {
    const std::size_t start = cid * chunk_size;
    if (start >= n) {
      continue;
    }
    const std::size_t end = std::min(n, start + chunk_size);
    out_interior.insert(
        out_interior.end(), base_interior.begin() + start, base_interior.begin() + end);
  }
}

void buildSnakeInteriorByInX(const std::vector<unsigned>& base_interior,
                             const std::vector<FlatScanChain::CellData>& cells,
                             std::vector<unsigned>& out_interior)
{
  out_interior = base_interior;
  if (out_interior.empty()) {
    return;
  }

  int min_x = cells[out_interior.front()].in_x;
  int max_x = min_x;
  for (const unsigned idx : out_interior) {
    min_x = std::min(min_x, cells[idx].in_x);
    max_x = std::max(max_x, cells[idx].in_x);
  }

  const double n = static_cast<double>(out_interior.size());
  const int target_bins = std::max(1, static_cast<int>(std::sqrt(n)));
  const int range = max_x - min_x;
  const int step = std::max(1, range / std::max(1, target_bins));

  std::sort(out_interior.begin(), out_interior.end(), [&](unsigned a, unsigned b) {
    const int ba = (cells[a].in_x - min_x) / step;
    const int bb = (cells[b].in_x - min_x) / step;
    if (ba != bb) {
      return ba < bb;
    }
    const bool flip = (ba & 1) != 0;
    if (!flip) {
      if (cells[a].in_y != cells[b].in_y) {
        return cells[a].in_y < cells[b].in_y;
      }
    } else {
      if (cells[a].in_y != cells[b].in_y) {
        return cells[a].in_y > cells[b].in_y;
      }
    }
    return a < b;
  });
}

void buildRandomProjectionInterior(const std::vector<unsigned>& base_interior,
                                   const std::vector<FlatScanChain::CellData>& cells,
                                   std::mt19937_64& rng,
                                   bool use_out,
                                   std::vector<unsigned>& out_interior)
{
  out_interior = base_interior;
  if (out_interior.empty()) {
    return;
  }
  std::uniform_int_distribution<int64_t> dist(-1024, 1024);
  int64_t a = 0;
  int64_t b = 0;
  for (int tries = 0; tries < 4; ++tries) {
    a = dist(rng);
    b = dist(rng);
    if (a != 0 || b != 0) {
      break;
    }
  }
  if (a == 0 && b == 0) {
    a = 1;
  }

  auto proj = [&](unsigned idx) -> int64_t {
    const auto& c = cells[idx];
    const int64_t x = use_out ? c.out_x : c.in_x;
    const int64_t y = use_out ? c.out_y : c.in_y;
    return a * x + b * y;
  };

  std::sort(out_interior.begin(), out_interior.end(), [&](unsigned x, unsigned y) {
    const int64_t px = proj(x);
    const int64_t py = proj(y);
    if (px != py) {
      return px < py;
    }
    return x < y;
  });
}

void buildStochasticGreedyPath(unsigned begin_idx,
                               unsigned end_idx,
                               const std::vector<unsigned>& base_interior,
                               const std::vector<FlatScanChain::CellData>& cells,
                               std::mt19937_64& rng,
                               std::size_t sample_size,
                               std::vector<unsigned>& out_path,
                               std::vector<unsigned>& scratch_remaining)
{
  out_path.clear();
  out_path.reserve(base_interior.size() + 2);
  out_path.push_back(begin_idx);

  scratch_remaining = base_interior;

  unsigned current = begin_idx;
  while (!scratch_remaining.empty()) {
    const std::size_t n = scratch_remaining.size();
    const std::size_t k = std::max<std::size_t>(1, std::min(sample_size, n));

    std::size_t best_pos = 0;
    int64_t best_cost = std::numeric_limits<int64_t>::max();
    for (std::size_t t = 0; t < k; ++t) {
      const std::size_t pos = static_cast<std::size_t>(rng() % n);
      const unsigned cand = scratch_remaining[pos];
      const int64_t cost = uclaEdgeCost(cells[current], cells[cand]);
      if (cost < best_cost
          || (cost == best_cost && cand < scratch_remaining[best_pos])) {
        best_cost = cost;
        best_pos = pos;
      }
    }

    const unsigned chosen = scratch_remaining[best_pos];
    scratch_remaining[best_pos] = scratch_remaining.back();
    scratch_remaining.pop_back();
    out_path.push_back(chosen);
    current = chosen;
  }

  out_path.push_back(end_idx);
}

void buildDeterministicGreedyPath(
    unsigned begin_idx,
    unsigned end_idx,
    const std::vector<unsigned>& base_interior,
    const std::vector<FlatScanChain::CellData>& cells,
    std::vector<unsigned>& out_path,
    std::vector<unsigned>& scratch_remaining)
{
  out_path.clear();
  out_path.reserve(base_interior.size() + 2);
  out_path.push_back(begin_idx);

  scratch_remaining = base_interior;

  unsigned current = begin_idx;
  while (!scratch_remaining.empty()) {
    std::size_t best_pos = 0;
    int64_t best_cost = std::numeric_limits<int64_t>::max();
    for (std::size_t pos = 0; pos < scratch_remaining.size(); ++pos) {
      const unsigned cand = scratch_remaining[pos];
      const int64_t cost = uclaEdgeCost(cells[current], cells[cand]);
      if (cost < best_cost
          || (cost == best_cost && cand < scratch_remaining[best_pos])) {
        best_cost = cost;
        best_pos = pos;
      }
    }

    const unsigned chosen = scratch_remaining[best_pos];
    scratch_remaining[best_pos] = scratch_remaining.back();
    scratch_remaining.pop_back();
    out_path.push_back(chosen);
    current = chosen;
  }

  out_path.push_back(end_idx);
}

struct AssignmentConstructorSpec
{
  int window = 12;
  std::size_t keep_per_row = 48;
  std::size_t break_candidates = 4;
  std::size_t fallback_per_column = 8;
};

struct AssignmentSortOrder
{
  std::vector<unsigned> nodes;
  std::vector<int> keys;
};

struct AssignmentArc
{
  unsigned left = 0;
  unsigned right = 0;
  int64_t cost = 0;
};

struct OpenedCycleOption
{
  std::vector<unsigned> path;
  unsigned tail = 0;
  unsigned head = 0;
  int64_t removed_cost = 0;
  int64_t retained_max = 0;
};

template <typename KeyFn>
AssignmentSortOrder buildAssignmentSortOrder(
    const std::vector<FlatScanChain::CellData>& cells,
    unsigned begin_idx,
    KeyFn key_fn)
{
  std::vector<std::pair<int, unsigned>> keyed_nodes;
  keyed_nodes.reserve(cells.size() - 1);
  for (unsigned node = 0; node < cells.size(); ++node) {
    if (node == begin_idx) {
      continue;
    }
    keyed_nodes.emplace_back(key_fn(node), node);
  }

  std::sort(keyed_nodes.begin(), keyed_nodes.end(), [](const auto& a, const auto& b) {
    if (a.first != b.first) {
      return a.first < b.first;
    }
    return a.second < b.second;
  });

  AssignmentSortOrder order;
  order.nodes.reserve(keyed_nodes.size());
  order.keys.reserve(keyed_nodes.size());
  for (const auto& [key, node] : keyed_nodes) {
    order.keys.push_back(key);
    order.nodes.push_back(node);
  }
  return order;
}

bool isLegalAssignmentArc(unsigned left,
                          unsigned right,
                          unsigned begin_idx,
                          unsigned end_idx,
                          std::size_t node_count)
{
  if (left >= node_count || right >= node_count) {
    return false;
  }
  if (left == right) {
    return false;
  }
  if (left == end_idx) {
    return right == begin_idx;
  }
  if (right == begin_idx) {
    return false;
  }
  if (left == begin_idx && right == end_idx && node_count > 2) {
    return false;
  }
  return true;
}

std::vector<std::vector<unsigned>> extractAssignmentCycles(
    const std::vector<unsigned>& successor)
{
  const std::size_t node_count = successor.size();
  std::vector<char> processed(node_count, 0);
  std::vector<int> local_pos(node_count, -1);
  std::vector<unsigned> walk;
  walk.reserve(node_count);
  std::vector<std::vector<unsigned>> cycles;

  for (unsigned start = 0; start < node_count; ++start) {
    if (processed[start]) {
      continue;
    }

    unsigned current = start;
    walk.clear();
    while (!processed[current] && local_pos[current] < 0) {
      local_pos[current] = static_cast<int>(walk.size());
      walk.push_back(current);
      current = successor[current];
    }

    if (local_pos[current] >= 0) {
      cycles.emplace_back(walk.begin() + local_pos[current], walk.end());
    }

    for (const unsigned node : walk) {
      processed[node] = 1;
      local_pos[node] = -1;
    }
  }

  return cycles;
}

int64_t cycleWorstEdgeCost(const std::vector<FlatScanChain::CellData>& cells,
                           const std::vector<unsigned>& cycle)
{
  int64_t worst = 0;
  for (std::size_t idx = 0; idx < cycle.size(); ++idx) {
    const unsigned left = cycle[idx];
    const unsigned right = cycle[(idx + 1) % cycle.size()];
    worst = std::max(worst, uclaEdgeCost(cells[left], cells[right]));
  }
  return worst;
}

std::vector<OpenedCycleOption> buildOpenedCycleOptions(
    const std::vector<FlatScanChain::CellData>& cells,
    const std::vector<unsigned>& cycle,
    std::size_t max_options)
{
  std::vector<std::pair<int64_t, std::size_t>> ranked_edges;
  ranked_edges.reserve(cycle.size());
  for (std::size_t edge_idx = 0; edge_idx < cycle.size(); ++edge_idx) {
    const unsigned tail = cycle[edge_idx];
    const unsigned head = cycle[(edge_idx + 1) % cycle.size()];
    ranked_edges.emplace_back(uclaEdgeCost(cells[tail], cells[head]), edge_idx);
  }

  std::sort(
      ranked_edges.begin(), ranked_edges.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
          return a.first > b.first;
        }
        return a.second < b.second;
      });

  const std::size_t keep = std::min(max_options, ranked_edges.size());
  std::vector<OpenedCycleOption> options;
  options.reserve(keep);
  for (std::size_t rank = 0; rank < keep; ++rank) {
    const std::size_t edge_idx = ranked_edges[rank].second;
    OpenedCycleOption option;
    option.tail = cycle[edge_idx];
    option.head = cycle[(edge_idx + 1) % cycle.size()];
    option.removed_cost = ranked_edges[rank].first;

    option.path.reserve(cycle.size());
    option.path.insert(
        option.path.end(), cycle.begin() + edge_idx + 1, cycle.end());
    option.path.insert(
        option.path.end(), cycle.begin(), cycle.begin() + edge_idx + 1);

    option.retained_max = 0;
    for (std::size_t pos = 0; pos + 1 < option.path.size(); ++pos) {
      option.retained_max = std::max(
          option.retained_max,
          uclaEdgeCost(cells[option.path[pos]], cells[option.path[pos + 1]]));
    }
    options.push_back(std::move(option));
  }

  return options;
}

bool isValidHamiltonianPath(const std::vector<unsigned>& path,
                            std::size_t node_count,
                            unsigned begin_idx,
                            unsigned end_idx)
{
  if (path.size() != node_count) {
    return false;
  }
  if (path.empty()) {
    return node_count == 0;
  }
  if (path.front() != begin_idx || path.back() != end_idx) {
    return false;
  }

  std::vector<char> seen(node_count, 0);
  for (const unsigned node : path) {
    if (node >= node_count || seen[node]) {
      return false;
    }
    seen[node] = 1;
  }
  return true;
}

bool buildAssignmentPatchedPath(const std::vector<FlatScanChain::CellData>& cells,
                                unsigned begin_idx,
                                unsigned end_idx,
                                std::vector<unsigned>& out_path)
{
  const std::size_t node_count = cells.size();
  if (node_count == 0) {
    out_path.clear();
    return true;
  }
  if (node_count <= 2) {
    out_path = {begin_idx, end_idx};
    return true;
  }

  const AssignmentSortOrder order_in_x = buildAssignmentSortOrder(
      cells, begin_idx, [&](unsigned node) { return cells[node].in_x; });
  const AssignmentSortOrder order_in_y = buildAssignmentSortOrder(
      cells, begin_idx, [&](unsigned node) { return cells[node].in_y; });
  const AssignmentSortOrder order_in_sum = buildAssignmentSortOrder(
      cells,
      begin_idx,
      [&](unsigned node) { return cells[node].in_x + cells[node].in_y; });
  const AssignmentSortOrder order_in_diff = buildAssignmentSortOrder(
      cells,
      begin_idx,
      [&](unsigned node) { return cells[node].in_x - cells[node].in_y; });

  constexpr std::array<AssignmentConstructorSpec, 4> specs{{
      {12, 48, 4, 8},
      {16, 64, 4, 8},
      {24, 96, 4, 8},
      {32, 128, 4, 8},
  }};

  std::vector<unsigned> incoming_count(node_count, 0);
  std::vector<uint32_t> candidate_marks(node_count, 0);
  std::vector<unsigned> row_candidates;
  row_candidates.reserve(160);
  std::vector<std::pair<int64_t, unsigned>> scored_candidates;
  scored_candidates.reserve(160);

  for (const auto& spec : specs) {
    std::vector<AssignmentArc> arcs;
    arcs.reserve((node_count - 1) * (spec.keep_per_row + 2));
    std::fill(incoming_count.begin(), incoming_count.end(), 0);

    arcs.push_back({end_idx, begin_idx, 0});
    incoming_count[begin_idx] = 1;

    uint32_t epoch = 0;
    bool row_failed = false;
    for (unsigned row = 0; row < node_count; ++row) {
      if (row == end_idx) {
        continue;
      }

      ++epoch;
      if (epoch == 0) {
        std::fill(candidate_marks.begin(), candidate_marks.end(), 0);
        ++epoch;
      }

      row_candidates.clear();
      auto add_from_order = [&](const AssignmentSortOrder& order, int query_key) {
        const auto it
            = std::lower_bound(order.keys.begin(), order.keys.end(), query_key);
        const std::size_t pos = static_cast<std::size_t>(
            std::distance(order.keys.begin(), it));
        const std::size_t radius = static_cast<std::size_t>(spec.window);
        const std::size_t lo = (pos > radius) ? pos - radius : 0;
        const std::size_t hi = std::min(order.nodes.size(), pos + radius + 1);
        for (std::size_t idx = lo; idx < hi; ++idx) {
          const unsigned column = order.nodes[idx];
          if (!isLegalAssignmentArc(
                  row, column, begin_idx, end_idx, node_count)) {
            continue;
          }
          if (candidate_marks[column] == epoch) {
            continue;
          }
          candidate_marks[column] = epoch;
          row_candidates.push_back(column);
        }
      };

      add_from_order(order_in_x, cells[row].out_x);
      add_from_order(order_in_y, cells[row].out_y);
      add_from_order(order_in_sum, cells[row].out_x + cells[row].out_y);
      add_from_order(order_in_diff, cells[row].out_x - cells[row].out_y);

      if (row != begin_idx && candidate_marks[end_idx] != epoch) {
        candidate_marks[end_idx] = epoch;
        row_candidates.push_back(end_idx);
      }

      scored_candidates.clear();
      for (const unsigned column : row_candidates) {
        scored_candidates.emplace_back(
            uclaEdgeCost(cells[row], cells[column]), column);
      }

      if (scored_candidates.empty()) {
        row_failed = true;
        break;
      }

      std::sort(scored_candidates.begin(),
                scored_candidates.end(),
                [](const auto& a, const auto& b) {
                  if (a.first != b.first) {
                    return a.first < b.first;
                  }
                  return a.second < b.second;
                });

      const std::size_t keep = std::min(spec.keep_per_row, scored_candidates.size());
      for (std::size_t idx = 0; idx < keep; ++idx) {
        const auto [cost, column] = scored_candidates[idx];
        arcs.push_back({row, column, cost});
        incoming_count[column] += 1;
      }
    }

    if (row_failed) {
      continue;
    }

    for (unsigned column = 0; column < node_count; ++column) {
      if (column == begin_idx || incoming_count[column] != 0) {
        continue;
      }

      std::vector<std::pair<int64_t, unsigned>> best_rows;
      best_rows.reserve(spec.fallback_per_column);
      for (unsigned row = 0; row < node_count; ++row) {
        if (!isLegalAssignmentArc(row, column, begin_idx, end_idx, node_count)) {
          continue;
        }
        best_rows.emplace_back(uclaEdgeCost(cells[row], cells[column]), row);
      }
      if (best_rows.empty()) {
        row_failed = true;
        break;
      }
      std::sort(best_rows.begin(), best_rows.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) {
          return a.first < b.first;
        }
        return a.second < b.second;
      });
      const std::size_t keep
          = std::min(spec.fallback_per_column, best_rows.size());
      for (std::size_t idx = 0; idx < keep; ++idx) {
        arcs.push_back({best_rows[idx].second, column, best_rows[idx].first});
      }
    }

    if (row_failed) {
      continue;
    }

    operations_research::SimpleLinearSumAssignment assignment;
    for (const auto& arc : arcs) {
      assignment.AddArcWithCost(
          static_cast<operations_research::NodeIndex>(arc.left),
          static_cast<operations_research::NodeIndex>(arc.right),
          static_cast<operations_research::CostValue>(arc.cost));
    }

    const auto status = assignment.Solve();
    if (status != operations_research::SimpleLinearSumAssignment::OPTIMAL) {
      continue;
    }

    std::vector<unsigned> successor(node_count, begin_idx);
    for (unsigned node = 0; node < node_count; ++node) {
      successor[node] = static_cast<unsigned>(assignment.RightMate(node));
    }
    if (successor[end_idx] != begin_idx) {
      continue;
    }

    std::vector<std::vector<unsigned>> cycles = extractAssignmentCycles(successor);
    std::size_t main_cycle_idx = cycles.size();
    for (std::size_t idx = 0; idx < cycles.size(); ++idx) {
      if (std::find(cycles[idx].begin(), cycles[idx].end(), begin_idx)
          != cycles[idx].end()) {
        main_cycle_idx = idx;
        break;
      }
    }
    if (main_cycle_idx == cycles.size()) {
      continue;
    }

    const std::vector<unsigned>& main_cycle = cycles[main_cycle_idx];
    const auto end_it = std::find(main_cycle.begin(), main_cycle.end(), end_idx);
    if (end_it == main_cycle.end()) {
      continue;
    }

    std::vector<unsigned> candidate_path;
    candidate_path.reserve(node_count);
    candidate_path.insert(candidate_path.end(), std::next(end_it), main_cycle.end());
    candidate_path.insert(candidate_path.end(), main_cycle.begin(), std::next(end_it));
    if (candidate_path.empty() || candidate_path.front() != begin_idx
        || candidate_path.back() != end_idx) {
      continue;
    }

    std::vector<std::vector<unsigned>> side_cycles;
    side_cycles.reserve(cycles.size() - 1);
    for (std::size_t idx = 0; idx < cycles.size(); ++idx) {
      if (idx == main_cycle_idx) {
        continue;
      }
      side_cycles.push_back(std::move(cycles[idx]));
    }
    std::sort(side_cycles.begin(),
              side_cycles.end(),
              [&](const auto& a, const auto& b) {
                const int64_t worst_a = cycleWorstEdgeCost(cells, a);
                const int64_t worst_b = cycleWorstEdgeCost(cells, b);
                if (worst_a != worst_b) {
                  return worst_a > worst_b;
                }
                return a.size() > b.size();
              });

    bool patch_failed = false;
    for (const auto& cycle : side_cycles) {
      const std::vector<OpenedCycleOption> options
          = buildOpenedCycleOptions(cells, cycle, spec.break_candidates);
      if (options.empty()) {
        patch_failed = true;
        break;
      }

      const OpenedCycleOption* best_option = nullptr;
      std::size_t best_insert_pos = 0;
      int64_t best_delta = std::numeric_limits<int64_t>::max();
      int64_t best_local_max = std::numeric_limits<int64_t>::max();
      int64_t best_removed_cost = -1;

      for (const auto& option : options) {
        for (std::size_t pos = 0; pos + 1 < candidate_path.size(); ++pos) {
          const unsigned left = candidate_path[pos];
          const unsigned right = candidate_path[pos + 1];
          const int64_t removed_path_cost = uclaEdgeCost(cells[left], cells[right]);
          const int64_t enter_cost = uclaEdgeCost(cells[left], cells[option.head]);
          const int64_t exit_cost = uclaEdgeCost(cells[option.tail], cells[right]);
          const int64_t delta = enter_cost + exit_cost - removed_path_cost;
          const int64_t local_max = std::max(
              option.retained_max, std::max(enter_cost, exit_cost));

          if (delta < best_delta
              || (delta == best_delta && local_max < best_local_max)
              || (delta == best_delta && local_max == best_local_max
                  && option.removed_cost > best_removed_cost)
              || (delta == best_delta && local_max == best_local_max
                  && option.removed_cost == best_removed_cost
                  && pos < best_insert_pos)) {
            best_delta = delta;
            best_local_max = local_max;
            best_removed_cost = option.removed_cost;
            best_insert_pos = pos;
            best_option = &option;
          }
        }
      }

      if (best_option == nullptr) {
        patch_failed = true;
        break;
      }

      candidate_path.insert(candidate_path.begin() + best_insert_pos + 1,
                            best_option->path.begin(),
                            best_option->path.end());
    }

    if (patch_failed) {
      continue;
    }
    if (!isValidHamiltonianPath(candidate_path, node_count, begin_idx, end_idx)) {
      continue;
    }

    out_path = std::move(candidate_path);
    return true;
  }

  out_path.clear();
  return false;
}

struct PortfolioRestartData
{
  std::vector<unsigned> interior_base;
  std::vector<unsigned> primary_path;
  std::vector<unsigned> greedy_path;
  std::array<std::vector<unsigned>, 4> deterministic_paths;
  std::vector<unsigned> assignment_path;
  std::vector<unsigned> best_chunk_base;
};

enum class PortfolioFamily : unsigned
{
  GreedyNearest = 0,
  RandomShuffle,
  SortedInX,
  SortedInY,
  SnakeInX,
  SortedOutX,
  ChunkShuffle,
  ProjectionIn,
  ProjectionOut,
  StochasticGreedy32,
  StochasticGreedy256,
  AssignmentCycleCover,
  Count
};

constexpr std::size_t kPortfolioFamilyCount
    = static_cast<std::size_t>(PortfolioFamily::Count);

struct PortfolioFamilyState
{
  uint64_t attempts = 0;
  bool has_score = false;
  UclaPathScore best_score;
};

std::size_t portfolioFamilyIndex(PortfolioFamily family)
{
  return static_cast<std::size_t>(family);
}

bool portfolioFamilyStateLess(const PortfolioFamilyState& a,
                              const PortfolioFamilyState& b)
{
  if (a.has_score != b.has_score) {
    return a.has_score && !b.has_score;
  }
  if (a.has_score) {
    if (uclaMetricLess(a.best_score, b.best_score)) {
      return true;
    }
    if (uclaMetricLess(b.best_score, a.best_score)) {
      return false;
    }
  }
  return a.attempts < b.attempts;
}

PortfolioFamily choosePortfolioFamily(
    const std::array<PortfolioFamilyState, kPortfolioFamilyCount>& family_states,
    uint64_t mix)
{
  std::array<std::size_t, kPortfolioFamilyCount> unexplored{};
  std::size_t unexplored_count = 0;
  std::array<std::size_t, kPortfolioFamilyCount> order{};
  for (std::size_t i = 0; i < kPortfolioFamilyCount; ++i) {
    order[i] = i;
    if (family_states[i].attempts == 0) {
      unexplored[unexplored_count++] = i;
    }
  }

  if (unexplored_count > 0) {
    const std::size_t pick
        = static_cast<std::size_t>(mix % static_cast<uint64_t>(unexplored_count));
    return static_cast<PortfolioFamily>(unexplored[pick]);
  }

  std::sort(order.begin(),
            order.end(),
            [&](std::size_t a, std::size_t b) {
              if (portfolioFamilyStateLess(family_states[a], family_states[b])) {
                return true;
              }
              if (portfolioFamilyStateLess(family_states[b], family_states[a])) {
                return false;
              }
              return a < b;
            });

  const uint64_t total_weight = static_cast<uint64_t>(kPortfolioFamilyCount)
                                * static_cast<uint64_t>(kPortfolioFamilyCount + 1)
                                / 2ULL;
  uint64_t ticket = mix % total_weight;
  for (std::size_t rank = 0; rank < kPortfolioFamilyCount; ++rank) {
    const uint64_t weight = static_cast<uint64_t>(kPortfolioFamilyCount - rank);
    if (ticket < weight) {
      return static_cast<PortfolioFamily>(order[rank]);
    }
    ticket -= weight;
  }

  return static_cast<PortfolioFamily>(order[0]);
}

PortfolioRestartData buildPortfolioRestartData(
    const std::vector<FlatScanChain::CellData>& cells,
    unsigned begin_idx,
    unsigned end_idx)
{
  PortfolioRestartData data;
  data.interior_base.reserve(cells.size() - 2);
  for (unsigned i = 1; i + 1 < cells.size(); ++i) {
    data.interior_base.push_back(i);
  }

  const std::vector<unsigned> sorted_in_x = preSortInterior(
      data.interior_base,
      [&](unsigned idx) {
        return std::make_pair(cells[idx].in_x, cells[idx].in_y);
      });
  const std::vector<unsigned> sorted_in_y = preSortInterior(
      data.interior_base,
      [&](unsigned idx) {
        return std::make_pair(cells[idx].in_y, cells[idx].in_x);
      });
  const std::vector<unsigned> sorted_out_x = preSortInterior(
      data.interior_base,
      [&](unsigned idx) {
        return std::make_pair(cells[idx].out_x, cells[idx].out_y);
      });

  std::vector<unsigned> restart_path;
  restart_path.reserve(cells.size());
  std::vector<unsigned> scratch_interior;
  scratch_interior.reserve(data.interior_base.size());
  std::vector<unsigned> scratch_remaining;
  scratch_remaining.reserve(data.interior_base.size());

  buildDeterministicGreedyPath(begin_idx,
                               end_idx,
                               data.interior_base,
                               cells,
                               restart_path,
                               scratch_remaining);
  data.greedy_path = restart_path;

  buildPathFromInterior(begin_idx, end_idx, sorted_in_x, restart_path);
  data.deterministic_paths[0] = restart_path;

  buildPathFromInterior(begin_idx, end_idx, sorted_in_y, restart_path);
  data.deterministic_paths[1] = restart_path;

  buildSnakeInteriorByInX(data.interior_base, cells, scratch_interior);
  buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
  data.deterministic_paths[2] = restart_path;

  buildPathFromInterior(begin_idx, end_idx, sorted_out_x, restart_path);
  data.deterministic_paths[3] = restart_path;

  buildAssignmentPatchedPath(cells, begin_idx, end_idx, data.assignment_path);
  if (!data.assignment_path.empty()) {
    uclaRefinePathL1(cells, data.assignment_path, /*end_fixed=*/true, std::nullopt);
    const auto assignment_polish_budget
        = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
            std::chrono::duration<double>(0.35));
    const auto assignment_polish_deadline
        = std::chrono::steady_clock::now() + assignment_polish_budget;
    uclaRefinePathL1Intensive(cells,
                              data.assignment_path,
                              /*end_fixed=*/true,
                              assignment_polish_deadline);
  }

  const std::vector<unsigned>* best_seed_path = &data.greedy_path;
  UclaPathScore best_seed_score = uclaEvaluatePath(cells, data.greedy_path);
  for (const auto& path : data.deterministic_paths) {
    const UclaPathScore score = uclaEvaluatePath(cells, path);
    if (uclaScoreLess(score, path, best_seed_score, *best_seed_path)) {
      best_seed_score = score;
      best_seed_path = &path;
    }
  }
  data.primary_path = *best_seed_path;
  if (data.primary_path.size() > 2) {
    data.best_chunk_base.assign(
        data.primary_path.begin() + 1, data.primary_path.end() - 1);
  }
  if (data.best_chunk_base.empty()) {
    data.best_chunk_base = data.interior_base;
  }

  return data;
}

bool buildPortfolioFamilyPath(
    PortfolioFamily family,
    unsigned begin_idx,
    unsigned end_idx,
    const PortfolioRestartData& restart_data,
    const std::vector<FlatScanChain::CellData>& cells,
    std::mt19937_64& rng,
    std::vector<unsigned>& restart_path,
    std::vector<unsigned>& scratch_interior,
    std::vector<unsigned>& scratch_remaining,
    bool& template_seed)
{
  template_seed = false;
  switch (family) {
    case PortfolioFamily::GreedyNearest:
      restart_path = restart_data.greedy_path;
      template_seed = !restart_path.empty();
      return template_seed;
    case PortfolioFamily::RandomShuffle:
      scratch_interior = restart_data.interior_base;
      std::shuffle(scratch_interior.begin(), scratch_interior.end(), rng);
      buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
      return true;
    case PortfolioFamily::SortedInX:
      restart_path = restart_data.deterministic_paths[0];
      template_seed = !restart_path.empty();
      return template_seed;
    case PortfolioFamily::SortedInY:
      restart_path = restart_data.deterministic_paths[1];
      template_seed = !restart_path.empty();
      return template_seed;
    case PortfolioFamily::SnakeInX:
      restart_path = restart_data.deterministic_paths[2];
      template_seed = !restart_path.empty();
      return template_seed;
    case PortfolioFamily::SortedOutX:
      restart_path = restart_data.deterministic_paths[3];
      template_seed = !restart_path.empty();
      return template_seed;
    case PortfolioFamily::ChunkShuffle: {
      const std::vector<unsigned>& chunk_base = restart_data.best_chunk_base.empty()
                                                    ? restart_data.interior_base
                                                    : restart_data.best_chunk_base;
      const std::size_t chunk = std::max<std::size_t>(32, chunk_base.size() / 64);
      buildChunkShuffledInterior(chunk_base, chunk, rng, scratch_interior);
      buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
      return true;
    }
    case PortfolioFamily::ProjectionIn:
      buildRandomProjectionInterior(
          restart_data.interior_base, cells, rng, false, scratch_interior);
      buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
      return true;
    case PortfolioFamily::ProjectionOut:
      buildRandomProjectionInterior(
          restart_data.interior_base, cells, rng, true, scratch_interior);
      buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
      return true;
    case PortfolioFamily::StochasticGreedy32:
      buildStochasticGreedyPath(begin_idx,
                                end_idx,
                                restart_data.interior_base,
                                cells,
                                rng,
                                32,
                                restart_path,
                                scratch_remaining);
      return true;
    case PortfolioFamily::StochasticGreedy256:
      buildStochasticGreedyPath(begin_idx,
                                end_idx,
                                restart_data.interior_base,
                                cells,
                                rng,
                                256,
                                restart_path,
                                scratch_remaining);
      return true;
    case PortfolioFamily::AssignmentCycleCover:
      if (restart_data.assignment_path.empty()) {
        return false;
      }
      restart_path = restart_data.assignment_path;
      template_seed = true;
      return true;
    case PortfolioFamily::Count:
      break;
  }
  return false;
}

void diversifyTemplatePath(std::vector<unsigned>& path,
                           uint64_t mix,
                           uint64_t prior_attempts,
                           std::mt19937_64& rng)
{
  if (path.size() < 8) {
    return;
  }
  if (prior_attempts == 0 && ((mix >> 22) & 1ULL) == 0ULL) {
    return;
  }

  const int bridge_kicks = 1 + static_cast<int>((mix >> 23) & 1ULL);
  for (int kick = 0; kick < bridge_kicks; ++kick) {
    doubleBridgeKick(path, rng);
  }
  if (((mix >> 24) & 1ULL) != 0ULL) {
    reverseSegmentKick(path, rng);
  }
}

MultiStartResult runMultiStartWorker(
    const std::vector<FlatScanChain::CellData>& cells,
    unsigned begin_idx,
    unsigned end_idx,
    const PortfolioRestartData& restart_data,
    const abkscanopt::Optimizer::Params& params,
    const UclaScanOptParams& ucla_params,
    const std::optional<std::chrono::steady_clock::time_point>& deadline,
    std::atomic<uint64_t>& next_restart,
    std::mutex& global_state_mu,
    std::atomic<int64_t>& global_best_cost_atomic,
    std::atomic<int64_t>& elite_worst_cost_atomic,
    int64_t& global_best_cost,
    std::vector<unsigned>& global_best_path,
    std::vector<EliteEntry>& global_elite,
    std::array<PortfolioFamilyState, kPortfolioFamilyCount>& family_states)
{
  MultiStartResult out;

  FlatScanChain chain(cells, begin_idx, end_idx);

  const std::vector<unsigned>& interior_base = restart_data.interior_base;
  std::vector<unsigned> restart_path;
  restart_path.reserve(cells.size());
  std::vector<unsigned> scratch_interior;
  scratch_interior.reserve(interior_base.size());
  std::vector<unsigned> scratch_remaining;
  scratch_remaining.reserve(interior_base.size());

  const unsigned base_major_loops = std::max(1U, params.majorLoops);
  const unsigned quick_major_loops = std::max(1U, base_major_loops / 10U);
  const unsigned medium_major_loops
      = std::max(quick_major_loops, (base_major_loops + 1U) / 2U);

  auto mulClamp = [](unsigned a, unsigned mul) -> unsigned {
    const uint64_t v = static_cast<uint64_t>(a) * mul;
    if (v > std::numeric_limits<unsigned>::max()) {
      return std::numeric_limits<unsigned>::max();
    }
    return static_cast<unsigned>(v);
  };
  const unsigned heavy_major_loops
      = std::max(base_major_loops, mulClamp(base_major_loops, 4));

  constexpr std::size_t kMaxElite = 8;

  while (true) {
    if (deadline.has_value()
        && std::chrono::steady_clock::now() >= deadline.value()) {
      break;
    }

    const uint64_t restart_id = next_restart.fetch_add(1);
    if (restart_id >= ucla_params.restarts) {
      break;
    }

    const uint64_t rand_seed64 = (restart_id == 0)
                                     ? ucla_params.seed
                                     : splitmix64(ucla_params.seed ^ restart_id);

    abkscanopt::Optimizer::Params local_params = params;
    PortfolioFamily selected_family = PortfolioFamily::GreedyNearest;
    uint64_t selected_family_prior_attempts = 0;
    bool selected_family_active = false;

    if (restart_id == 0) {
      restart_path = restart_data.primary_path.empty() ? restart_data.greedy_path
                                                       : restart_data.primary_path;
      if (restart_path.empty()) {
        buildDeterministicGreedyPath(begin_idx,
                                     end_idx,
                                     interior_base,
                                     cells,
                                     restart_path,
                                     scratch_remaining);
      }
      chain.restorePath(restart_path);
      local_params.majorLoops = base_major_loops;
    } else {
      std::mt19937_64 rng(splitmix64(rand_seed64));
      const uint64_t mix = splitmix64(rand_seed64);

      const bool have_elite
          = elite_worst_cost_atomic.load(std::memory_order_relaxed)
                != std::numeric_limits<int64_t>::max()
            || global_best_cost_atomic.load(std::memory_order_relaxed)
                   != std::numeric_limits<int64_t>::max();
      const bool super_intensify = have_elite && ((mix & 63ULL) == 0ULL);
      const bool intensify = have_elite && ((mix & 15ULL) < 4ULL);

      if (intensify) {
        std::vector<unsigned> parent_a;
        std::vector<unsigned> parent_b;
        std::vector<unsigned> start_path;
        {
          std::lock_guard<std::mutex> lock(global_state_mu);
          if (global_elite.size() >= 2 && ((mix >> 7) & 1ULL) != 0ULL) {
            const std::size_t k = std::min<std::size_t>(4, global_elite.size());
            const std::size_t pick_a = static_cast<std::size_t>(rng() % k);
            std::size_t pick_b = static_cast<std::size_t>(rng() % k);
            if (pick_b == pick_a && k > 1) {
              pick_b = (pick_b + 1) % k;
            }
            parent_a = global_elite[pick_a].path;
            parent_b = global_elite[pick_b].path;
          } else if (!global_elite.empty()) {
            const std::size_t pick = static_cast<std::size_t>(
                rng() % std::min<std::size_t>(4, global_elite.size()));
            start_path = global_elite[pick].path;
          } else {
            start_path = global_best_path;
          }
        }

        if (!parent_a.empty() && !parent_b.empty()) {
          if (!orderCrossoverOX(parent_a, parent_b, rng, start_path)) {
            start_path = std::move(parent_a);
          }
        }
        if (start_path.empty()) {
          start_path = restart_data.primary_path;
        }

        if (!start_path.empty()) {
          const int kicks = (super_intensify ? 2 : 1)
                            + static_cast<int>((mix >> 4) & 3ULL);
          for (int kick = 0; kick < kicks; ++kick) {
            doubleBridgeKick(start_path, rng);
          }
          if (super_intensify || (((mix >> 6) & 1ULL) != 0ULL)) {
            reverseSegmentKick(start_path, rng);
          }
          chain.restorePath(start_path);
          local_params.majorLoops = super_intensify ? heavy_major_loops
                                                    : base_major_loops;
          if (super_intensify) {
            local_params.nnear = std::min(128U, local_params.nnear * 2U);
            local_params.nDescents = std::min(16U, local_params.nDescents * 2U);
          }
        } else {
          local_params.majorLoops = quick_major_loops;
        }
      } else {
        bool template_seed = false;
        constexpr std::size_t kAdaptiveFamilyMinInterior = 2500;
        if (restart_data.interior_base.size() >= kAdaptiveFamilyMinInterior) {
          std::lock_guard<std::mutex> lock(global_state_mu);
          selected_family = choosePortfolioFamily(family_states, mix);
          auto& family_state = family_states[portfolioFamilyIndex(selected_family)];
          selected_family_prior_attempts = family_state.attempts;
          family_state.attempts += 1;
          selected_family_active = true;
        } else {
          const bool use_assignment_mode
              = !restart_data.assignment_path.empty()
                && (((mix >> 28) & 15ULL) == 0ULL);
          if (use_assignment_mode) {
            selected_family = PortfolioFamily::AssignmentCycleCover;
          } else {
            const unsigned mode = static_cast<unsigned>(mix % 9ULL);
            switch (mode) {
              case 0:
                selected_family = PortfolioFamily::RandomShuffle;
                break;
              case 1:
                selected_family = PortfolioFamily::SortedInX;
                break;
              case 2:
                selected_family = PortfolioFamily::SortedInY;
                break;
              case 3:
                selected_family = PortfolioFamily::SnakeInX;
                break;
              case 4:
                selected_family = PortfolioFamily::GreedyNearest;
                break;
              case 5:
                selected_family = PortfolioFamily::ChunkShuffle;
                break;
              case 6:
                selected_family = (((mix >> 16) & 1ULL) != 0ULL)
                                      ? PortfolioFamily::ProjectionOut
                                      : PortfolioFamily::ProjectionIn;
                break;
              case 7:
                selected_family = (((mix >> 20) & 1ULL) != 0ULL)
                                      ? PortfolioFamily::StochasticGreedy256
                                      : PortfolioFamily::StochasticGreedy32;
                break;
              case 8:
              default:
                selected_family = PortfolioFamily::SortedOutX;
                break;
            }
          }
        }

        if (buildPortfolioFamilyPath(selected_family,
                                     begin_idx,
                                     end_idx,
                                     restart_data,
                                     cells,
                                     rng,
                                     restart_path,
                                     scratch_interior,
                                     scratch_remaining,
                                     template_seed)) {
          if (template_seed
              && (selected_family == PortfolioFamily::AssignmentCycleCover
                  || ((mix >> 21) & 1ULL) != 0ULL)) {
            diversifyTemplatePath(
                restart_path, mix, selected_family_prior_attempts, rng);
          }
        } else {
          scratch_interior = interior_base;
          std::shuffle(scratch_interior.begin(), scratch_interior.end(), rng);
          buildPathFromInterior(begin_idx, end_idx, scratch_interior, restart_path);
          selected_family = PortfolioFamily::RandomShuffle;
          template_seed = false;
        }

        chain.restorePath(restart_path);
        local_params.majorLoops = quick_major_loops;
        if (selected_family == PortfolioFamily::AssignmentCycleCover) {
          local_params.majorLoops = medium_major_loops;
        }

        if (((mix >> 25) & 1ULL) != 0ULL) {
          local_params.only2Opt = true;
          local_params.nnear = std::max(8U, local_params.nnear / 2U);
          local_params.nDescents = std::max(1U, local_params.nDescents / 2U);
        }
      }
    }

    RandomRawUnsigned randuns(clampToDeterministicSeed(rand_seed64));
    optimizeChainSafe(chain, randuns, local_params, deadline);

    if (!uclaTimeExpired(deadline)) {
      restart_path = chain.getPath();
      uclaRefinePathL1(cells, restart_path, /*end_fixed=*/true, deadline);
      chain.restorePath(restart_path);
    }

    const UclaPathScore final_score = uclaEvaluatePath(cells, chain.getPath());
    const int64_t cost = final_score.cost;
    if (cost < out.best_cost
        || (cost == out.best_cost && lexLess(chain.getPath(), out.best_path))) {
      out.best_cost = cost;
      out.best_path = chain.getPath();
    }

    {
      std::lock_guard<std::mutex> lock(global_state_mu);
      if (selected_family_active) {
        auto& family_state = family_states[portfolioFamilyIndex(selected_family)];
        if (!family_state.has_score
            || uclaMetricLess(final_score, family_state.best_score)) {
          family_state.has_score = true;
          family_state.best_score = final_score;
        }
      }

      if (global_best_path.empty() || cost < global_best_cost
          || (cost == global_best_cost
              && lexLess(chain.getPath(), global_best_path))) {
        global_best_cost = cost;
        global_best_path = chain.getPath();
        global_best_cost_atomic.store(cost, std::memory_order_relaxed);
      }

      updateElitePool(global_elite, kMaxElite, cost, chain.getPath());
      if (global_elite.empty() || global_elite.size() < kMaxElite) {
        elite_worst_cost_atomic.store(std::numeric_limits<int64_t>::max(),
                                      std::memory_order_relaxed);
      } else {
        elite_worst_cost_atomic.store(global_elite.back().cost,
                                      std::memory_order_relaxed);
      }
    }
  }

  return out;
}


}  // namespace

std::vector<std::size_t> UclaScanOptOrder(
    const std::vector<std::string_view>& names,
    const std::vector<std::pair<int, int>>& scan_in_pts,
    const std::vector<std::pair<int, int>>& scan_out_pts,
    const std::pair<int, int>& begin,
    const std::pair<int, int>& end,
    const UclaScanOptParams& params)
{
  const std::size_t n = names.size();
  if (scan_in_pts.size() != n || scan_out_pts.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: size mismatch");
  }
  if (n <= 1) {
    return n == 1 ? std::vector<std::size_t>{0} : std::vector<std::size_t>{};
  }

  std::vector<FlatScanChain::CellData> cells;
  cells.reserve(n + 2);

  {
    FlatScanChain::CellData c;
    c.in_x = begin.first;
    c.in_y = begin.second;
    c.out_x = begin.first;
    c.out_y = begin.second;
    c.name = "__begin__";
    cells.push_back(std::move(c));
  }

  for (std::size_t i = 0; i < n; ++i) {
    FlatScanChain::CellData c;
    c.in_x = scan_in_pts[i].first;
    c.in_y = scan_in_pts[i].second;
    c.out_x = scan_out_pts[i].first;
    c.out_y = scan_out_pts[i].second;
    c.name = std::string(names[i]);
    cells.push_back(std::move(c));
  }

  {
    FlatScanChain::CellData c;
    c.in_x = end.first;
    c.in_y = end.second;
    c.out_x = end.first;
    c.out_y = end.second;
    c.name = "__end__";
    cells.push_back(std::move(c));
  }

  const unsigned begin_idx = 0;
  const unsigned end_idx = static_cast<unsigned>(cells.size() - 1);

  FlatScanChain chain(cells, begin_idx, end_idx);

  configureSeedHandlerOnce(params.seed);

  abkscanopt::Optimizer::Params p;
  p.majorLoops = clampToUnsigned(params.major_loops);
  p.nDescents = clampToUnsigned(params.n_descents);
  p.kickMove = clampToUnsigned(params.kick_move);
  p.nnear = clampToUnsigned(params.n_near);
  p.only2Opt = params.only_2opt;
  p.temp_control = params.temp_control;
  p.zeroTemp = params.zero_temp;

  RandomRawUnsigned randuns(clampToDeterministicSeed(params.seed));
  optimizeChainSafe(chain, randuns, p, std::nullopt);

  const auto& path = chain.getPath();
  if (path.size() != n + 2) {
    throw std::runtime_error("UclaScanOptOrder: unexpected path length");
  }
  if (path.front() != begin_idx || path.back() != end_idx) {
    throw std::runtime_error("UclaScanOptOrder: endpoints moved unexpectedly");
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  for (std::size_t k = 1; k + 1 < path.size(); ++k) {
    const unsigned idx = path[k];
    if (idx == begin_idx || idx == end_idx) {
      throw std::runtime_error("UclaScanOptOrder: internal endpoint");
    }
    order.push_back(static_cast<std::size_t>(idx - 1));
  }

  if (order.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: ordering size mismatch");
  }

  return order;
}

std::vector<std::size_t> UclaScanOptOrderPortfolio(
    const std::vector<std::string_view>& names,
    const std::vector<std::pair<int, int>>& scan_in_pts,
    const std::vector<std::pair<int, int>>& scan_out_pts,
    const std::pair<int, int>& begin,
    const std::pair<int, int>& end,
    const UclaScanOptParams& params)
{
  const std::size_t n = names.size();
  if (scan_in_pts.size() != n || scan_out_pts.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: size mismatch");
  }
  if (n <= 1) {
    return n == 1 ? std::vector<std::size_t>{0} : std::vector<std::size_t>{};
  }

  std::vector<FlatScanChain::CellData> cells;
  cells.reserve(n + 2);

  {
    FlatScanChain::CellData c;
    c.in_x = begin.first;
    c.in_y = begin.second;
    c.out_x = begin.first;
    c.out_y = begin.second;
    c.name = "__begin__";
    cells.push_back(std::move(c));
  }

  for (std::size_t i = 0; i < n; ++i) {
    FlatScanChain::CellData c;
    c.in_x = scan_in_pts[i].first;
    c.in_y = scan_in_pts[i].second;
    c.out_x = scan_out_pts[i].first;
    c.out_y = scan_out_pts[i].second;
    c.name = std::string(names[i]);
    cells.push_back(std::move(c));
  }

  {
    FlatScanChain::CellData c;
    c.in_x = end.first;
    c.in_y = end.second;
    c.out_x = end.first;
    c.out_y = end.second;
    c.name = "__end__";
    cells.push_back(std::move(c));
  }

  const unsigned begin_idx = 0;
  const unsigned end_idx = static_cast<unsigned>(cells.size() - 1);

  FlatScanChain chain(cells, begin_idx, end_idx);

  configureSeedHandlerOnce(params.seed);

  abkscanopt::Optimizer::Params p;
  p.majorLoops = clampToUnsigned(params.major_loops);
  p.nDescents = clampToUnsigned(params.n_descents);
  p.kickMove = clampToUnsigned(params.kick_move);
  p.nnear = clampToUnsigned(params.n_near);
  p.only2Opt = params.only_2opt;
  p.temp_control = params.temp_control;
  p.zeroTemp = params.zero_temp;

  UclaScanOptParams effective = params;
  if (params.time_limit_seconds > 0.0 && effective.restarts <= 1) {
    effective.restarts = std::numeric_limits<uint64_t>::max();
  }

  const unsigned threads = std::max(1U, clampToUnsigned(effective.threads));
  const uint64_t restarts64 = std::max<uint64_t>(1, effective.restarts);

  const bool use_multistart = (restarts64 > 1);
  const bool use_threads = (threads > 1);
  const bool use_time_limit = (effective.time_limit_seconds > 0.0);

  if (!use_multistart && !use_time_limit && !use_threads) {
    RandomRawUnsigned randuns(clampToDeterministicSeed(params.seed));
    optimizeChainSafe(chain, randuns, p, std::nullopt);
  } else {
    std::optional<std::chrono::steady_clock::time_point> deadline;
    if (use_time_limit) {
      const auto dur
          = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
              std::chrono::duration<double>(effective.time_limit_seconds));
      deadline = std::chrono::steady_clock::now() + dur;
    }

    std::optional<std::chrono::steady_clock::time_point> search_deadline = deadline;
    if (deadline.has_value()) {
      const double total_s = std::max(0.0, effective.time_limit_seconds);
      const double reserve_s = std::clamp(0.1 * total_s, 1.0, 60.0);
      const auto reserve
          = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
              std::chrono::duration<double>(reserve_s));
      const auto now = std::chrono::steady_clock::now();
      if (now + reserve < deadline.value()) {
        search_deadline = deadline.value() - reserve;
      }
    }

    const unsigned workers
        = static_cast<unsigned>(std::min<uint64_t>(threads, restarts64));
    const PortfolioRestartData restart_data
        = buildPortfolioRestartData(cells, begin_idx, end_idx);
    std::vector<MultiStartResult> results;
    results.reserve(workers);

    std::atomic<uint64_t> next_restart{0};
    std::mutex global_state_mu;
    std::atomic<int64_t> global_best_cost_atomic{
        std::numeric_limits<int64_t>::max()};
    std::atomic<int64_t> elite_worst_cost_atomic{
        std::numeric_limits<int64_t>::max()};
    int64_t global_best_cost = std::numeric_limits<int64_t>::max();
    std::vector<unsigned> global_best_path;
    global_best_path.reserve(cells.size());
    std::vector<EliteEntry> global_elite;
    global_elite.reserve(8);
    std::array<PortfolioFamilyState, kPortfolioFamilyCount> family_states;

    auto seed_global_state = [&](const std::vector<unsigned>& path) {
      if (path.empty()) {
        return;
      }
      if (!isValidHamiltonianPath(path, cells.size(), begin_idx, end_idx)) {
        return;
      }
      if (std::any_of(global_elite.begin(),
                      global_elite.end(),
                      [&](const EliteEntry& entry) { return entry.path == path; })) {
        return;
      }

      const int64_t seed_cost = uclaPathCost(cells, path);
      if (global_best_path.empty() || seed_cost < global_best_cost
          || (seed_cost == global_best_cost && lexLess(path, global_best_path))) {
        global_best_cost = seed_cost;
        global_best_path = path;
      }
      updateElitePool(global_elite, 8, seed_cost, path);
    };

    seed_global_state(restart_data.primary_path);
    seed_global_state(restart_data.greedy_path);
    for (const auto& path : restart_data.deterministic_paths) {
      seed_global_state(path);
    }
    if (!global_best_path.empty()) {
      global_best_cost_atomic.store(global_best_cost, std::memory_order_relaxed);
    }
    if (global_elite.empty() || global_elite.size() < 8) {
      elite_worst_cost_atomic.store(std::numeric_limits<int64_t>::max(),
                                    std::memory_order_relaxed);
    } else {
      elite_worst_cost_atomic.store(global_elite.back().cost,
                                    std::memory_order_relaxed);
    }
    std::vector<std::thread> pool;
    pool.reserve(workers);

    for (unsigned i = 0; i < workers; ++i) {
      results.emplace_back();
      pool.emplace_back([&, i]() {
        results[i] = runMultiStartWorker(
            cells,
            begin_idx,
            end_idx,
            restart_data,
            p,
            effective,
            search_deadline,
            next_restart,
            global_state_mu,
            global_best_cost_atomic,
            elite_worst_cost_atomic,
            global_best_cost,
            global_best_path,
            global_elite,
            family_states);
      });
    }
    for (auto& t : pool) {
      t.join();
    }

    MultiStartResult best;
    for (const auto& r : results) {
      if (r.best_path.empty()) {
        continue;
      }
      if (r.best_cost < best.best_cost
          || (r.best_cost == best.best_cost
              && std::lexicographical_compare(r.best_path.begin(),
                                              r.best_path.end(),
                                              best.best_path.begin(),
                                              best.best_path.end()))) {
        best = r;
      }
    }

    if (!best.best_path.empty()) {
      std::vector<unsigned> refined = best.best_path;

      if (deadline.has_value() && !uclaTimeExpired(deadline)) {
        std::vector<EliteEntry> elite_copy;
        {
          std::lock_guard<std::mutex> lock(global_state_mu);
          elite_copy = global_elite;
        }
        if (elite_copy.empty()) {
          elite_copy.push_back(EliteEntry{uclaPathCost(cells, refined), refined});
        }

        std::optional<std::chrono::steady_clock::time_point> elite_deadline
            = deadline;
        if (deadline.has_value()) {
          const auto now = std::chrono::steady_clock::now();
          const double remaining_s = std::max(
              0.0, std::chrono::duration<double>(deadline.value() - now).count());
          const double exact_reserve_s = std::clamp(0.33 * remaining_s, 0.5, 2.0);
          const auto exact_reserve = std::chrono::duration_cast<
              std::chrono::steady_clock::duration>(
              std::chrono::duration<double>(exact_reserve_s));
          if (now + exact_reserve < deadline.value()) {
            elite_deadline = deadline.value() - exact_reserve;
          }
        }

        const unsigned refine_workers = workers;
        std::vector<std::thread> refine_pool;
        refine_pool.reserve(refine_workers);
        std::vector<EliteEntry> refined_out;
        refined_out.resize(refine_workers);
        for (unsigned t = 0; t < refine_workers; ++t) {
          refine_pool.emplace_back([&, t]() {
            std::mt19937_64 rng(
                splitmix64(effective.seed ^ (0x9e3779b97f4a7c15ULL * (t + 1))));
            EliteEntry best_local;
            best_local.path = elite_copy[t % elite_copy.size()].path;
            best_local.cost = uclaPathCost(cells, best_local.path);

            while (!uclaTimeExpired(elite_deadline)) {
              std::vector<unsigned> cand = best_local.path;
              const uint64_t mix = splitmix64(rng());
              const int kicks = 1 + static_cast<int>(mix & 3ULL);
              for (int k = 0; k < kicks; ++k) {
                doubleBridgeKick(cand, rng);
              }
              if (((mix >> 2) & 1ULL) != 0ULL) {
                reverseSegmentKick(cand, rng);
              }

              const auto slice = std::chrono::duration_cast<
                  std::chrono::steady_clock::duration>(
                  std::chrono::duration<double>(0.25));
              const auto round_deadline = std::min(
                  elite_deadline.value_or(std::chrono::steady_clock::now() + slice),
                  std::chrono::steady_clock::now() + slice);

              uclaRefinePathL1Intensive(
                  cells, cand, /*end_fixed=*/true, round_deadline);
              const int64_t cost = uclaPathCost(cells, cand);
              if (cost < best_local.cost
                  || (cost == best_local.cost
                      && lexLess(cand, best_local.path))) {
                best_local.cost = cost;
                best_local.path.swap(cand);
              }
            }

            refined_out[t] = std::move(best_local);
          });
        }
        for (auto& t : refine_pool) {
          t.join();
        }

        EliteEntry best_elite;
        for (const auto& e : refined_out) {
          if (e.path.empty()) {
            continue;
          }
          if (e.cost < best_elite.cost
              || (e.cost == best_elite.cost
                  && lexLess(e.path, best_elite.path))) {
            best_elite = e;
          }
        }
        if (!best_elite.path.empty()) {
          refined = std::move(best_elite.path);
          if (!uclaTimeExpired(deadline)) {
            std::vector<EliteEntry> exact_out;
            exact_out.resize(refine_workers);
            std::vector<std::thread> exact_pool;
            exact_pool.reserve(refine_workers);
            for (unsigned t = 0; t < refine_workers; ++t) {
              exact_pool.emplace_back([&, t]() {
                std::mt19937_64 rng(splitmix64(
                    effective.seed ^ (0xd1b54a32d192ed03ULL + t)));
                std::vector<unsigned> cand = refined;
                UclaPathScore best_score = uclaEvaluatePath(cells, cand);
                std::vector<unsigned> best_path = cand;

                for (int pass = 0; pass < 8; ++pass) {
                  if (uclaTimeExpired(deadline)) {
                    break;
                  }
                  const bool improved = uclaRefinePathExactWindows(
                      cells, cand, /*end_fixed=*/true, rng, deadline);
                  if (!uclaTimeExpired(deadline)) {
                    uclaRefinePathL1Intensive(
                        cells, cand, /*end_fixed=*/true, deadline);
                  }
                  if (!uclaTimeExpired(deadline)) {
                    uclaRefinePathJumpPlateau(
                        cells, cand, /*end_fixed=*/true, deadline);
                  }

                  const UclaPathScore cand_score = uclaEvaluatePath(cells, cand);
                  if (uclaScoreLess(cand_score, cand, best_score, best_path)) {
                    best_score = cand_score;
                    best_path = cand;
                  }
                  if (!improved) {
                    break;
                  }
                }

                exact_out[t] = EliteEntry{uclaPathCost(cells, best_path), best_path};
              });
            }
            for (auto& t : exact_pool) {
              t.join();
            }

            UclaPathScore refined_score = uclaEvaluatePath(cells, refined);
            for (const auto& e : exact_out) {
              if (e.path.empty()) {
                continue;
              }
              const UclaPathScore cand_score = uclaEvaluatePath(cells, e.path);
              if (uclaScoreLess(cand_score, e.path, refined_score, refined)) {
                refined = e.path;
                refined_score = cand_score;
              }
            }
            if (!uclaTimeExpired(deadline)) {
              uclaRefinePathJumpPlateau(cells, refined, /*end_fixed=*/true, deadline);
            }
          }
        }
      } else {
        const auto tiny
            = std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(0.10));
        const auto polish_deadline = std::chrono::steady_clock::now() + tiny;
        std::mt19937_64 polish_rng(splitmix64(effective.seed ^ 0x51f15e5dULL));
        uclaRefinePathL1Intensive(
            cells, refined, /*end_fixed=*/true, polish_deadline);
        if (!uclaTimeExpired(polish_deadline)) {
          uclaRefinePathExactWindows(
              cells, refined, /*end_fixed=*/true, polish_rng, polish_deadline);
        }
        if (!uclaTimeExpired(polish_deadline)) {
          uclaRefinePathJumpPlateau(
              cells, refined, /*end_fixed=*/true, polish_deadline);
        }
      }

      chain.restorePath(refined);
    }
  }

  const auto& path = chain.getPath();
  if (path.size() != n + 2) {
    throw std::runtime_error("UclaScanOptOrder: unexpected path length");
  }
  if (path.front() != begin_idx || path.back() != end_idx) {
    throw std::runtime_error("UclaScanOptOrder: endpoints moved unexpectedly");
  }

  std::vector<std::size_t> order;
  order.reserve(n);
  for (std::size_t k = 1; k + 1 < path.size(); ++k) {
    const unsigned idx = path[k];
    if (idx == begin_idx || idx == end_idx) {
      throw std::runtime_error("UclaScanOptOrder: internal endpoint");
    }
    order.push_back(static_cast<std::size_t>(idx - 1));
  }

  if (order.size() != n) {
    throw std::runtime_error("UclaScanOptOrder: ordering size mismatch");
  }

  return order;
}

}  // namespace dft
