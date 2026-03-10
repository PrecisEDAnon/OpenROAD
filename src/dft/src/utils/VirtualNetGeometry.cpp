// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "VirtualNetGeometry.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "odb/dbShape.h"

namespace dft {
namespace {

constexpr int64_t kInfDistance = std::numeric_limits<int64_t>::max() / 8;

int64_t manhattanDist(const odb::Point& a,
                      const odb::Point& b,
                      double vertical_weight)
{
  const int64_t dx = std::abs(static_cast<int64_t>(a.x())
                              - static_cast<int64_t>(b.x()));
  const int64_t dy = std::abs(static_cast<int64_t>(a.y())
                              - static_cast<int64_t>(b.y()));
  const int64_t wy = static_cast<int64_t>(std::llround(vertical_weight * dy));
  return dx + wy;
}

int64_t manhattanPointToRectDist(const odb::Point& p,
                                 const odb::Rect& r,
                                 double vertical_weight)
{
  int64_t dx = 0;
  if (p.x() < r.xMin()) {
    dx = static_cast<int64_t>(r.xMin()) - p.x();
  } else if (p.x() > r.xMax()) {
    dx = static_cast<int64_t>(p.x()) - r.xMax();
  }

  int64_t dy = 0;
  if (p.y() < r.yMin()) {
    dy = static_cast<int64_t>(r.yMin()) - p.y();
  } else if (p.y() > r.yMax()) {
    dy = static_cast<int64_t>(p.y()) - r.yMax();
  }
  const int64_t wy = static_cast<int64_t>(std::llround(vertical_weight * dy));
  return dx + wy;
}

odb::Point closestPointOnRect(const odb::Point& p, const odb::Rect& r)
{
  const int x = std::clamp(p.x(), r.xMin(), r.xMax());
  const int y = std::clamp(p.y(), r.yMin(), r.yMax());
  return odb::Point(x, y);
}

bool betterVirtualPinCandidate(int64_t dist_a,
                               const odb::Point& pt_a,
                               VirtualPinKind kind_a,
                               int64_t dist_b,
                               const odb::Point& pt_b,
                               VirtualPinKind kind_b)
{
  if (dist_a != dist_b) {
    return dist_a < dist_b;
  }
  if (pt_a.x() != pt_b.x()) {
    return pt_a.x() < pt_b.x();
  }
  if (pt_a.y() != pt_b.y()) {
    return pt_a.y() < pt_b.y();
  }
  // Prefer terminals over boxes as a deterministic tie-break.
  return static_cast<int>(kind_a) > static_cast<int>(kind_b);
}

}  // namespace

VirtualSubnet buildVirtualSubnetFromNet(odb::dbNet* net)
{
  VirtualSubnet geom;
  if (net == nullptr) {
    return geom;
  }

  for (odb::dbGuide* guide : net->getGuides()) {
    geom.boxes.push_back(guide->getBox());
  }

  if (geom.boxes.empty()) {
    if (odb::dbWire* wire = net->getWire()) {
      odb::dbWireShapeItr itr;
      odb::dbShape shape;
      for (itr.begin(wire); itr.next(shape);) {
        geom.boxes.push_back(shape.getBox());
      }
    }
  }

  if (geom.boxes.empty()) {
    for (odb::dbITerm* iterm : net->getITerms()) {
      const odb::Rect bbox = iterm->getBBox();
      geom.terminals.emplace_back(bbox.xMin(), bbox.yMin());
    }
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const odb::Rect bbox = bterm->getBBox();
      geom.terminals.emplace_back(bbox.xMin(), bbox.yMin());
    }
  }

  return geom;
}

VirtualSubnet clipVirtualSubnet(const VirtualSubnet& subnet, const odb::Rect& clip)
{
  VirtualSubnet out;

  for (const odb::Rect& box : subnet.boxes) {
    const int x0 = std::max(box.xMin(), clip.xMin());
    const int y0 = std::max(box.yMin(), clip.yMin());
    const int x1 = std::min(box.xMax(), clip.xMax());
    const int y1 = std::min(box.yMax(), clip.yMax());
    if (x0 <= x1 && y0 <= y1) {
      out.boxes.emplace_back(x0, y0, x1, y1);
    }
  }

  for (const odb::Point& p : subnet.terminals) {
    if (p.x() >= clip.xMin() && p.x() <= clip.xMax() && p.y() >= clip.yMin()
        && p.y() <= clip.yMax()) {
      out.terminals.push_back(p);
    }
  }

  return out;
}

VirtualPinResult distanceToVirtualSubnet(const odb::Point& pin,
                                         const VirtualSubnet& subnet,
                                         double vertical_weight)
{
  VirtualPinResult best;
  best.dist = kInfDistance;
  best.closest = odb::Point(0, 0);
  best.kind = VirtualPinKind::None;

  for (const odb::Rect& box : subnet.boxes) {
    const int64_t d = manhattanPointToRectDist(pin, box, vertical_weight);
    const odb::Point closest = closestPointOnRect(pin, box);
    if (!best.valid()
        || betterVirtualPinCandidate(
            d, closest, VirtualPinKind::Box, best.dist, best.closest, best.kind)) {
      best.dist = d;
      best.closest = closest;
      best.kind = VirtualPinKind::Box;
      if (best.dist == 0) {
        // We can’t do better than 0, but keep scanning in case another 0
        // produces a deterministically smaller closest point.
      }
    }
  }

  for (const odb::Point& term : subnet.terminals) {
    const int64_t d = manhattanDist(pin, term, vertical_weight);
    if (!best.valid()
        || betterVirtualPinCandidate(d,
                                     term,
                                     VirtualPinKind::Terminal,
                                     best.dist,
                                     best.closest,
                                     best.kind)) {
      best.dist = d;
      best.closest = term;
      best.kind = VirtualPinKind::Terminal;
      if (best.dist == 0) {
        // Keep scanning for a smaller terminal point, if any.
      }
    }
  }

  return best;
}

}  // namespace dft

