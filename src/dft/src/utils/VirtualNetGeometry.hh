// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#pragma once

#include <cstdint>
#include <vector>

#include "odb/db.h"
#include "odb/geom.h"

namespace dft {

// A virtual subnet is a lightweight geometric proxy for a (possibly-routed)
// net, represented as a union of axis-aligned rectangles plus terminal points.
// It can be sourced from an OpenDB net's guides/wire-shapes, or constructed
// synthetically (future work) for scan/timing/power optimization.
struct VirtualSubnet
{
  std::vector<odb::Rect> boxes;
  std::vector<odb::Point> terminals;

  bool empty() const { return boxes.empty() && terminals.empty(); }
};

enum class VirtualPinKind
{
  None,
  Box,
  Terminal
};

struct VirtualPinResult
{
  int64_t dist = 0;
  odb::Point closest{0, 0};
  VirtualPinKind kind = VirtualPinKind::None;

  bool valid() const { return kind != VirtualPinKind::None; }
};

// Build a subnet geometry proxy from an OpenDB net:
// - prefers route guides, else route wire shapes, else pin bbox lower-left
//   corners.
VirtualSubnet buildVirtualSubnetFromNet(odb::dbNet* net);

// Clip a subnet geometry to a rectangular region (intersection for boxes,
// containment for terminals). Useful for constructing local "subnets" around a
// chosen virtual pin.
VirtualSubnet clipVirtualSubnet(const VirtualSubnet& subnet, const odb::Rect& clip);

// Compute the minimum weighted Manhattan distance from a point to a subnet and
// return the corresponding closest point (a "virtual pin") on that subnet.
VirtualPinResult distanceToVirtualSubnet(const odb::Point& pin,
                                         const VirtualSubnet& subnet,
                                         double vertical_weight);

}  // namespace dft

