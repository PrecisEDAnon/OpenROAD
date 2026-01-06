// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include "timingBase.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "db_sta/dbNetwork.hh"
#include "odb/db.h"
#include "nesterovBase.h"
#include "placerBase.h"
#include "rsz/Resizer.hh"
#include "sta/Fuzzy.hh"
#include "sta/Graph.hh"
#include "utl/Logger.h"

namespace gpl {

using utl::GPL;

// TimingBase
TimingBase::TimingBase() = default;

TimingBase::TimingBase(std::shared_ptr<NesterovBaseCommon> nbc,
                       rsz::Resizer* rs,
                       utl::Logger* log)
    : TimingBase()
{
  rs_ = rs;
  nbc_ = std::move(nbc);
  log_ = log;
}

void TimingBase::initTimingOverflowChk()
{
  timingOverflowChk_.clear();
  timingOverflowChk_.resize(timingNetWeightOverflow_.size(), false);
}

bool TimingBase::isTimingNetWeightOverflow(float overflow)
{
  int intOverflow = std::round(overflow * 100);
  // exception case handling
  if (timingNetWeightOverflow_.empty()
      || intOverflow > timingNetWeightOverflow_[0]) {
    return false;
  }

  bool needTdRun = false;
  for (int i = 0; i < timingNetWeightOverflow_.size(); i++) {
    if (timingNetWeightOverflow_[i] > intOverflow) {
      if (!timingOverflowChk_[i]) {
        timingOverflowChk_[i] = true;
        needTdRun = true;
      }
      continue;
    }
    return needTdRun;
  }
  return needTdRun;
}

void TimingBase::addTimingNetWeightOverflow(int overflow)
{
  std::vector<int>::iterator it = std::find(timingNetWeightOverflow_.begin(),
                                            timingNetWeightOverflow_.end(),
                                            overflow);

  // only push overflow when the overflow is not in vector.
  if (it == timingNetWeightOverflow_.end()) {
    timingNetWeightOverflow_.push_back(overflow);
  }

  // do sort in reverse order
  std::sort(timingNetWeightOverflow_.begin(),
            timingNetWeightOverflow_.end(),
            std::greater<int>());
}

void TimingBase::setTimingNetWeightOverflows(std::vector<int>& overflows)
{
  // sort by decreasing order
  std::sort(overflows.begin(), overflows.end(), std::greater<int>());
  for (auto& overflow : overflows) {
    addTimingNetWeightOverflow(overflow);
  }
  initTimingOverflowChk();
}

void TimingBase::deleteTimingNetWeightOverflow(int overflow)
{
  std::vector<int>::iterator it = std::find(timingNetWeightOverflow_.begin(),
                                            timingNetWeightOverflow_.end(),
                                            overflow);
  // only erase overflow when the overflow is in vector.
  if (it != timingNetWeightOverflow_.end()) {
    timingNetWeightOverflow_.erase(it);
  }
}

void TimingBase::clearTimingNetWeightOverflow()
{
  timingNetWeightOverflow_.clear();
}

size_t TimingBase::getTimingNetWeightOverflowSize() const
{
  return timingNetWeightOverflow_.size();
}

void TimingBase::setTimingNetWeightMax(float max)
{
  net_weight_max_ = max;
}

void TimingBase::setDynamicWeightOptions(const DynamicWeightOptions& options)
{
  dynamic_options_ = options;
}

void TimingBase::setCorridorOptions(const CorridorOptions& options)
{
  corridor_options_ = options;
}

void TimingBase::setSpringOptions(const SpringOptions& options)
{
  spring_options_ = options;
}

void TimingBase::setNesterovBases(const std::vector<std::shared_ptr<NesterovBase>>& bases)
{
  nbs_ = bases;
  clearCorridorMask();
  clearPathSprings();
  locator_dirty_ = true;
}

void TimingBase::setECPScale(float scale)
{
  ecp_scale_ = scale;
}

void TimingBase::setECPWeightMax(float max)
{
  ecp_weight_max_ = max;
}

void TimingBase::setECPTopEndpointFrac(float frac)
{
  ecp_top_endpoint_frac_ = frac;
}

void TimingBase::updateResizerGuards(const NesterovPlaceVars& vars)
{
  if (!rs_) return;
  rs_->setNonCriticalSlackBudget(vars.td_noncrit_slack_budget_ns);
  rs_->setSetupGuardCapNs(vars.td_setup_guard_cap_ns);
  rs_->setSetupGuardWindowNs(vars.td_guard_window_ns);
  rs_->setSetupSlackGuard(vars.td_setup_guard_cap_ns);
  rs_->setPostCtsHoldFloorNs(vars.td_post_cts_hold_floor_ns);
  rs_->setRebufferCloneGateFanout(vars.td_rebuffer_clone_gate_fanout);
  rs_->setCloneGroupFanout(vars.td_clone_group_fanout);
  rs_->setGrPickRadiusTiles(vars.td_gr_pick_radius);
  guard_defaults_.noncrit_budget = vars.td_noncrit_slack_budget_ns;
  guard_defaults_.setup_guard_cap = vars.td_setup_guard_cap_ns;
  guard_defaults_.guard_window = vars.td_guard_window_ns;
  guard_defaults_.setup_guard_floor = vars.td_setup_guard_cap_ns;
  guard_defaults_.valid = true;
  applyResizerGuards(SlackRegime::LowSlack);
}

void TimingBase::rebuildGCellLocator()
{
  gcell_locator_.clear();
  for (size_t nb_idx = 0; nb_idx < nbs_.size(); ++nb_idx) {
    auto nb = nbs_[nb_idx];
    const auto& handles = nb->getGCells();
    for (size_t i = 0; i < handles.size(); ++i) {
      const GCell* gCell = static_cast<const GCell*>(handles[i]);
      gcell_locator_[gCell] = CellLocator{nb_idx, i};
    }
  }
  locator_dirty_ = false;
}

std::optional<TimingBase::CellLocator> TimingBase::locateGCell(const GCell* gCell)
{
  if (locator_dirty_) rebuildGCellLocator();
  auto it = gcell_locator_.find(gCell);
  if (it != gcell_locator_.end()) return it->second;
  return std::nullopt;
}

TimingBase::SlackRegime TimingBase::classifyRegime(float wns)
{
  if (wns >= 0.05e-9) return SlackRegime::HighSlack; // 50ps
  if (wns >= -0.05e-9) return SlackRegime::MediumSlack;
  return SlackRegime::LowSlack;
}

std::vector<TimingBase::PathChain> TimingBase::collectPathChains(int top_count, float endpoint_slack_threshold, int min_length)
{
  std::vector<PathChain> chains;

  if (!rs_ || top_count <= 0) {
    return chains;
  }

  if (min_length <= 0) {
    min_length = 1;
  }
  
  auto staPaths = rs_->findWorstPaths(top_count, endpoint_slack_threshold);
  
  for (auto& staPath : staPaths) {
    PathChain chain;
    for (auto* vertex : staPath) {
      // Skip port vertices if they don't map to instances well?
      // We can check if the pin is a port of the block/top cell?
      // Or just try to map to instance.
      
      sta::Pin* pin = vertex->pin();
      if (!pin) continue;
      
      odb::dbITerm* iterm;
      odb::dbBTerm* bterm;
      odb::dbModITerm* moditerm;
      rs_->getDbNetwork()->staToDb(pin, iterm, bterm, moditerm);
      
      GCell* gCell = nullptr;
      if (iterm && iterm->getInst()) {
        gCell = nbc_->dbToNb(iterm->getInst());
      }
      
      if (gCell) {
        // Avoid duplicates in chain (e.g. input and output pins of same instance)
        if (chain.empty() || chain.back() != gCell) {
          chain.push_back(gCell);
        }
      }
    }
    
    if (chain.size() >= min_length) {
      chains.push_back(chain);
    }
  }
  
  return chains;
}

void TimingBase::applyCorridorMask(const std::vector<PathChain>& chains)
{
  for (auto& nb : nbs_) {
    int binCntX = nb->getBinCntX();
    int binCntY = nb->getBinCntY();
    int min_x = nb->getBinGrid().lx();
    int min_y = nb->getBinGrid().ly();
    double binSizeX = nb->getBinSizeX();
    double binSizeY = nb->getBinSizeY();
    
    std::vector<float> charge(binCntX * binCntY, 0.0f);
    
    for (const auto& chain : chains) {
      if (chain.size() < 2) continue;
      
      for (size_t i = 0; i < chain.size() - 1; ++i) {
        GCell* u = chain[i];
        GCell* v = chain[i+1];
        
        int x1 = u->cx();
        int y1 = u->cy();
        int x2 = v->cx();
        int y2 = v->cy();
        
        // Simple line walking
        float dist = std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
        int steps = std::max(1, (int)(dist / std::min(binSizeX, binSizeY)));
        
        for (int s = 0; s <= steps; ++s) {
          float t = (float)s / steps;
          int cx = x1 + t * (x2 - x1);
          int cy = y1 + t * (y2 - y1);
          
          int bx = (cx - min_x) / binSizeX;
          int by = (cy - min_y) / binSizeY;
          
          if (bx >= 0 && bx < binCntX && by >= 0 && by < binCntY) {
            // Apply to neighborhood
            int w = corridor_options_.width_bins;
            for (int dx = -w; dx <= w; ++dx) {
              for (int dy = -w; dy <= w; ++dy) {
                int nx = bx + dx;
                int ny = by + dy;
                if (nx >= 0 && nx < binCntX && ny >= 0 && ny < binCntY) {
                  charge[ny * binCntX + nx] += 1.0f;
                }
              }
            }
          }
        }
      }
    }
    
    std::vector<float> mask(binCntX * binCntY);
    for (size_t i = 0; i < charge.size(); ++i) {
       mask[i] = 1.0f / (1.0f + corridor_options_.charge_k * charge[i]);
    }
    nb->setCorridorMask(mask);
  }
}

void TimingBase::clearCorridorMask()
{
  for (auto& nb : nbs_) {
    nb->clearCorridorMask();
  }
}

void TimingBase::applyPathSprings(const std::vector<PathChain>& chains)
{
  std::vector<std::vector<FloatPoint>> all_forces(nbs_.size());
  for (size_t i = 0; i < nbs_.size(); ++i) {
    all_forces[i].resize(nbs_[i]->getGCells().size(), FloatPoint(0, 0));
  }
  
  for (const auto& chain : chains) {
    if (chain.size() < 2) continue;
    
    for (size_t i = 0; i < chain.size() - 1; ++i) {
      GCell* u = chain[i];
      GCell* v = chain[i+1];
      
      int dx = v->cx() - u->cx();
      int dy = v->cy() - u->cy();
      float dist = std::sqrt(dx*dx + dy*dy);
      
      if (dist < 1e-6) continue;
      
      float fx = spring_options_.k * dx / dist;
      float fy = spring_options_.k * dy / dist;
      
      if (auto loc_u = locateGCell(u)) {
         all_forces[loc_u->nb_index][loc_u->cell_index].x += fx;
         all_forces[loc_u->nb_index][loc_u->cell_index].y += fy;
      }
      if (auto loc_v = locateGCell(v)) {
         all_forces[loc_v->nb_index][loc_v->cell_index].x -= fx;
         all_forces[loc_v->nb_index][loc_v->cell_index].y -= fy;
      }
    }
  }
  
  for (size_t i = 0; i < nbs_.size(); ++i) {
    nbs_[i]->setPathSpringForces(all_forces[i]);
  }
}

void TimingBase::clearPathSprings()
{
  for (auto& nb : nbs_) {
    nb->clearPathSpringForces();
  }
}

bool TimingBase::crossesCongestedBins(GNet* gNet)
{
  if (nbs_.empty()) return false;
  
  // Simplified: check first region.
  auto nb = nbs_[0];
  auto& bg = nb->getBinGrid();
  auto& bins = nb->getBins();
  int binCntX = nb->getBinCntX();
  int binCntY = nb->getBinCntY();
  
  int lx = gNet->lx();
  int ly = gNet->ly();
  int ux = gNet->ux();
  int uy = gNet->uy();
  
  int min_bx = (lx - bg.lx()) / bg.getBinSizeX();
  int min_by = (ly - bg.ly()) / bg.getBinSizeY();
  int max_bx = (ux - bg.lx()) / bg.getBinSizeX();
  int max_by = (uy - bg.ly()) / bg.getBinSizeY();
  
  min_bx = std::max(0, std::min(binCntX - 1, min_bx));
  min_by = std::max(0, std::min(binCntY - 1, min_by));
  max_bx = std::max(0, std::min(binCntX - 1, max_bx));
  max_by = std::max(0, std::min(binCntY - 1, max_by));
  
  int visited = 0;
  int hot = 0;
  float total_overflow = 0;
  
  for (int bx = min_bx; bx <= max_bx; ++bx) {
    for (int by = min_by; by <= max_by; ++by) {
      int idx = by * binCntX + bx;
      float ov = bins[idx].getDensity() - bins[idx].getTargetDensity();
      if (ov > 0) {
        total_overflow += ov;
        if (ov > dynamic_options_.hot_bin_threshold) hot++;
      }
      visited++;
    }
  }
  
  if (visited == 0) return false;
  
  float avg_overflow = total_overflow / visited;
  float hot_frac = (float)hot / visited;
  
  return (avg_overflow >= dynamic_options_.congestion_gate) || 
         (hot_frac > dynamic_options_.hot_bin_fraction);
}

void TimingBase::applyECPWeights(SlackRegime regime)
{
  if (rs_->countECP(ecp_scale_) == 0) return;
  
  std::unordered_map<const sta::Net*, float> net2crit;
  rs_->collectECPNets(ecp_scale_, ecp_top_endpoint_frac_, net2crit);
  
  for (auto& gNet : nbc_->getGNets()) {
    sta::Net* staNet = gNet->getPbNet()->getDbNet() ? rs_->getDbNetwork()->dbToSta(gNet->getPbNet()->getDbNet()) : nullptr;
    if (!staNet) continue;
    
    if (net2crit.find(staNet) == net2crit.end()) continue;
    
    float crit = net2crit[staNet];
    float norm = std::min(1.0f, crit / 10.0f); // Normalize
    
    float base_w = 1.0f + (ecp_weight_max_ - 1.0f) * norm;
    
    if (regime == SlackRegime::HighSlack) {
       gNet->setTimingWeight(base_w);
    } else if (regime == SlackRegime::MediumSlack) {
       float old = gNet->getTimingWeight();
       float blended = 0.6f * old + 0.4f * base_w;
       gNet->setTimingWeight(blended);
    }
    // LowSlack: ignore ECP, focus on WNS
  }
}

void TimingBase::resetAllNetWeights()
{
  if (!nbc_) {
    return;
  }
  for (auto& gNet : nbc_->getGNets()) {
    gNet->setTimingWeight(1.0f);
  }
}

void TimingBase::applyResizerGuards(SlackRegime regime)
{
  if (!rs_ || !guard_defaults_.valid) {
    return;
  }

  auto clamp_non_neg = [](double value) { return value < 0.0 ? 0.0 : value; };

  double noncrit_budget = guard_defaults_.noncrit_budget;
  double guard_cap = guard_defaults_.setup_guard_cap;
  double guard_window = guard_defaults_.guard_window;
  double guard_floor = guard_defaults_.setup_guard_floor;

  switch (regime) {
    case SlackRegime::LowSlack:
      break;
    case SlackRegime::MediumSlack:
      noncrit_budget *= 0.5;
      guard_cap *= 0.8;
      guard_window *= 1.2;
      guard_floor *= 0.8;
      break;
    case SlackRegime::HighSlack:
      noncrit_budget = 0.0;
      guard_cap *= 0.5;
      guard_window *= 2.0;
      guard_floor *= 0.5;
      break;
  }

  rs_->setNonCriticalSlackBudget(clamp_non_neg(noncrit_budget));
  rs_->setSetupGuardCapNs(clamp_non_neg(guard_cap));
  rs_->setSetupGuardWindowNs(clamp_non_neg(guard_window));
  rs_->setSetupSlackGuard(clamp_non_neg(guard_floor));
}

bool TimingBase::executeTimingDriven(int iter,
                                     int /*total_iter*/,
                                     float avg_overflow,
                                     bool run_journal_restore)
{
  if (!rs_ || !nbc_) {
    return false;
  }

  // In virtual (journal-restore) iterations we only need timing analysis for
  // guidance; skip the expensive/perturbing repair_design step and avoid
  // growing the design when we won't keep changes anyway.
  rs_->findResizeSlacks(run_journal_restore, /*repair_design=*/!run_journal_restore);

  if (!run_journal_restore) {
    nbc_->fixPointers();
  }

  sta::NetSeq worst_slack_nets = rs_->resizeWorstSlackNets();
  if (worst_slack_nets.empty()) {
    log_->warn(GPL, 700, "Timing-driven: no net slacks found.");
    return false;
  }

  auto wns_opt = rs_->resizeNetSlack(worst_slack_nets.front());
  if (!wns_opt || std::isinf(wns_opt.value())) {
    log_->warn(GPL, 701, "Timing-driven: invalid worst slack measurement.");
    return false;
  }
  const float wns = wns_opt.value();

  SlackRegime regime = classifyRegime(wns);
  const bool dynamic_enabled = dynamic_options_.enable;
  const bool dynamic_update_due
      = dynamic_enabled
        && (dynamic_options_.update_period <= 1
            || ((iter + 1) % dynamic_options_.update_period) == 0);

  bool weights_updated = false;
  const float endpoint_threshold = std::isfinite(dynamic_options_.endpoint_slack_threshold)
                                       ? dynamic_options_.endpoint_slack_threshold
                                       : std::numeric_limits<float>::infinity();

  if (dynamic_update_due) {
    auto apply_dynamic = [&]() -> bool {
      sta::dbNetwork* db_network = rs_->getDbNetwork();
      if (db_network == nullptr) {
        return false;
      }

      std::vector<odb::dbNet*> candidate_db_nets;
      candidate_db_nets.reserve(worst_slack_nets.size());
      std::unordered_set<odb::dbNet*> seen_db_nets;
      auto add_candidate = [&](odb::dbNet* db_net) {
        if (db_net && seen_db_nets.insert(db_net).second) {
          candidate_db_nets.push_back(db_net);
        }
      };

      if (dynamic_options_.top_endpoints > 0) {
        auto critical_nets = rs_->criticalPathNets(dynamic_options_.top_endpoints,
                                                   endpoint_threshold);
        for (odb::dbNet* db_net : critical_nets) {
          add_candidate(db_net);
        }
      }

      for (const sta::Net* sta_net : worst_slack_nets) {
        odb::dbNet* db_net = db_network->staToDb(sta_net);
        add_candidate(db_net);
      }

      if (candidate_db_nets.empty()) {
        return false;
      }

      std::vector<std::pair<odb::dbNet*, float>> ordered;
      ordered.reserve(candidate_db_nets.size());
      std::unordered_map<odb::dbNet*, float> slack_map;
      for (odb::dbNet* db_net : candidate_db_nets) {
        if (auto slack_opt = rs_->resizeNetSlack(db_net)) {
          const float slack_value = slack_opt.value();
          if (std::isinf(slack_value)) {
            continue;
          }
          ordered.emplace_back(db_net, slack_value);
          slack_map[db_net] = slack_value;
        }
      }

      if (ordered.empty()) {
        return false;
      }

      std::sort(ordered.begin(),
                ordered.end(),
                [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

      const float slack_min = ordered.front().second;
      const float slack_max = ordered.back().second;

      float severity_ratio = 0.0f;
      if (dynamic_options_.severity_slack_norm > 0.0) {
        severity_ratio
            = std::clamp(static_cast<float>(-slack_min / dynamic_options_.severity_slack_norm),
                         0.0f,
                         dynamic_options_.severity_ratio_cap);
      }

      float ramp = dynamic_options_.ramp_iterations > 0
                       ? std::min(1.0f,
                                  static_cast<float>(iter + 1)
                                      / static_cast<float>(dynamic_options_.ramp_iterations))
                       : 1.0f;
      if (regime == SlackRegime::MediumSlack) {
        ramp *= 0.6f;
      }

      float w_max_eff = dynamic_options_.weight_max
                        * (1.0f + dynamic_options_.severity_weight_scale * severity_ratio);
      w_max_eff = std::min(dynamic_options_.severity_weight_limit, w_max_eff);
      w_max_eff = std::max(dynamic_options_.weight_min, w_max_eff);
      const float dynamic_span = std::max(0.0f, w_max_eff - dynamic_options_.weight_min);

      float coverage_base = dynamic_options_.initial_nets_percent
                            + (dynamic_options_.final_nets_percent
                               - dynamic_options_.initial_nets_percent)
                                  * ramp;
      float coverage_boost
          = std::min(dynamic_options_.severity_coverage_limit,
                     1.0f + dynamic_options_.severity_coverage_scale * severity_ratio);
      float coverage_fraction
          = std::clamp(coverage_base * coverage_boost / 100.0f, 0.0f, 1.0f);
      const double strong_target
          = std::ceil(static_cast<double>(ordered.size())
                      * static_cast<double>(coverage_fraction));
      size_t strong_count = std::max<size_t>(1, static_cast<size_t>(strong_target));

      std::unordered_set<odb::dbNet*> strong_nets;
      for (size_t i = 0; i < ordered.size() && i < strong_count; ++i) {
        strong_nets.insert(ordered[i].first);
      }

      float overflow_factor = 1.0f;
      if (dynamic_options_.overflow_limit > 0.0f) {
        const float normalized
            = std::clamp(avg_overflow / dynamic_options_.overflow_limit, 0.0f, 1.0f);
        overflow_factor
            = std::max(0.0f, 1.0f - dynamic_options_.congestion_alpha * normalized);
      }

      resetAllNetWeights();

      bool assigned = false;
      for (auto* gNet : nbc_->getGNets()) {
        if (!gNet || gNet->getGPins().size() <= 1) {
          continue;
        }

        Net* pb_net = gNet->getPbNet();
        if (!pb_net) {
          continue;
        }
        odb::dbNet* db_net = pb_net->getDbNet();
        if (!db_net) {
          continue;
        }

        auto slack_it = slack_map.find(db_net);
        if (slack_it == slack_map.end()) {
          continue;
        }

        if (crossesCongestedBins(gNet)) {
          continue;
        }

        const float slack = slack_it->second;
        const float slack_severity
            = dynamic_options_.slack_norm > 0.0
                  ? std::clamp(
                        static_cast<float>(-slack / dynamic_options_.slack_norm), 0.0f, 1.0f)
                  : 1.0f;

        const float hpwl
            = static_cast<float>((gNet->ux() - gNet->lx()) + (gNet->uy() - gNet->ly()));
        const float length_factor = dynamic_options_.length_norm > 0.0
                                        ? std::clamp(
                                              hpwl / static_cast<float>(dynamic_options_.length_norm),
                                              0.0f,
                                              1.0f)
                                        : 1.0f;

        const bool strong = strong_nets.find(db_net) != strong_nets.end();
        float weight = 1.0f;
        if (strong) {
          weight = dynamic_options_.weight_min
                   + dynamic_span * ramp * slack_severity * length_factor * overflow_factor;
          weight = std::clamp(weight, dynamic_options_.weight_min, w_max_eff);
        } else {
          const float denom = std::max(1e-9f, slack_max - slack_min);
          const float interp = std::clamp((slack_max - slack) / denom, 0.0f, 1.0f);
          const float fallback_max = std::min(w_max_eff, 1.0f + 0.25f * dynamic_span);
          weight = 1.0f + (fallback_max - 1.0f) * interp * ramp;
          weight = std::clamp(weight, 1.0f, fallback_max);
        }

        gNet->setTimingWeight(weight);
        assigned = true;
      }

      return assigned;
    };

    weights_updated = apply_dynamic();
  }

  if (!dynamic_options_.enable) {
    auto apply_legacy = [&]() -> bool {
      const auto slack_min_opt = rs_->resizeNetSlack(worst_slack_nets.front());
      const auto slack_max_opt = rs_->resizeNetSlack(worst_slack_nets.back());
      if (!slack_min_opt || !slack_max_opt) {
        return false;
      }

      const float slack_min = slack_min_opt.value();
      const float slack_max = slack_max_opt.value();
      if (std::isinf(slack_min) || std::isinf(slack_max)) {
        return false;
      }

      resetAllNetWeights();
      const float denom = std::max(1e-9f, slack_max - slack_min);
      int weighted = 0;
      for (auto* gNet : nbc_->getGNets()) {
        if (!gNet || gNet->getGPins().size() <= 1) {
          continue;
        }

        Net* pb_net = gNet->getPbNet();
        if (!pb_net) {
          continue;
        }
        odb::dbNet* db_net = pb_net->getDbNet();
        if (!db_net) {
          continue;
        }

        auto slack_opt = rs_->resizeNetSlack(db_net);
        if (!slack_opt) {
          continue;
        }
        const float slack = slack_opt.value();
        if (std::isinf(slack) || slack >= slack_max) {
          continue;
        }

        float weight = 1.0f;
        if (slack_max != slack_min) {
          weight = 1.0f
                   + (net_weight_max_ - 1.0f) * (slack_max - slack) / denom;
        }
        gNet->setTimingWeight(weight);
        ++weighted;
      }
      return weighted > 0;
    };

    weights_updated = apply_legacy();
  }

  if (weights_updated) {
    weights_initialized_ = true;
  } else if (!weights_initialized_) {
    resetAllNetWeights();
  }

  if (regime != SlackRegime::LowSlack) {
    applyECPWeights(regime);
  }

  const auto apply_spatial_guidance = [&](SlackRegime active_regime) {
    const float threshold = endpoint_threshold;

    if (active_regime == SlackRegime::LowSlack) {
      if (corridor_options_.enable) {
        auto chains = collectPathChains(corridor_options_.top_endpoints,
                                        threshold,
                                        corridor_options_.min_path_length);
        applyCorridorMask(chains);
      } else {
        clearCorridorMask();
      }

      if (spring_options_.enable && avg_overflow < spring_options_.overflow_gate) {
        auto chains = collectPathChains(spring_options_.top_paths,
                                        threshold,
                                        spring_options_.min_path_length);
        applyPathSprings(chains);
      } else {
        clearPathSprings();
      }
    } else if (active_regime == SlackRegime::MediumSlack) {
      clearCorridorMask();
      if (spring_options_.enable && avg_overflow < 0.1f) {
        auto chains = collectPathChains(spring_options_.top_paths,
                                        threshold,
                                        spring_options_.min_path_length);
        applyPathSprings(chains);
      } else {
        clearPathSprings();
      }
    } else {
      clearCorridorMask();
      clearPathSprings();
    }
  };

  apply_spatial_guidance(regime);
  applyResizerGuards(regime);

  return true;
}

// Keep the old signature for compatibility or if called from somewhere else, but redirect
bool TimingBase::executeTimingDriven(bool run_journal_restore) {
    return executeTimingDriven(0, 1000, 0.1f, run_journal_restore);
}

}  // namespace gpl
