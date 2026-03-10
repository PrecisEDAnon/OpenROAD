// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, The OpenROAD Authors

#include "VirtualNetGeometry.hh"

#include "gtest/gtest.h"

namespace dft {
namespace {

TEST(VirtualNetGeometry, DistanceAndClosestPoint)
{
  VirtualSubnet subnet;
  subnet.boxes.emplace_back(10, 10, 20, 20);
  subnet.terminals.emplace_back(0, 0);

  {
    const odb::Point p(15, 12);
    const VirtualPinResult res = distanceToVirtualSubnet(p, subnet, 1.0);
    EXPECT_TRUE(res.valid());
    EXPECT_EQ(res.dist, 0);
    EXPECT_EQ(res.closest.x(), p.x());
    EXPECT_EQ(res.closest.y(), p.y());
  }

  {
    const odb::Point p(5, 12);
    const VirtualPinResult res = distanceToVirtualSubnet(p, subnet, 1.0);
    EXPECT_TRUE(res.valid());
    EXPECT_EQ(res.dist, 5);
    EXPECT_EQ(res.closest.x(), 10);
    EXPECT_EQ(res.closest.y(), 12);
  }

  // Terminal beats the box when it is closer.
  {
    const odb::Point p(5, 0);
    const VirtualPinResult res = distanceToVirtualSubnet(p, subnet, 1.0);
    EXPECT_TRUE(res.valid());
    EXPECT_EQ(res.dist, 5);
    EXPECT_EQ(res.closest.x(), 0);
    EXPECT_EQ(res.closest.y(), 0);
    EXPECT_EQ(res.kind, VirtualPinKind::Terminal);
  }
}

TEST(VirtualNetGeometry, ClipSubnet)
{
  VirtualSubnet subnet;
  subnet.boxes.emplace_back(10, 10, 20, 20);
  subnet.terminals.emplace_back(0, 0);
  subnet.terminals.emplace_back(18, 19);

  const odb::Rect clip(15, 15, 25, 25);
  const VirtualSubnet clipped = clipVirtualSubnet(subnet, clip);

  ASSERT_EQ(clipped.boxes.size(), 1U);
  EXPECT_EQ(clipped.boxes[0].xMin(), 15);
  EXPECT_EQ(clipped.boxes[0].yMin(), 15);
  EXPECT_EQ(clipped.boxes[0].xMax(), 20);
  EXPECT_EQ(clipped.boxes[0].yMax(), 20);

  ASSERT_EQ(clipped.terminals.size(), 1U);
  EXPECT_EQ(clipped.terminals[0].x(), 18);
  EXPECT_EQ(clipped.terminals[0].y(), 19);
}

}  // namespace
}  // namespace dft

