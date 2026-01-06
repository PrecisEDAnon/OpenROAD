// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "rsz/Resizer.hh"

namespace rsz {
class Resizer;
}

namespace utl {
class Logger;
}

namespace gpl {

class NesterovBaseCommon;
class NesterovBase;
class NesterovPlaceVars;
class GNet;
class GCell;

class TimingBase
{
 public:
  struct DynamicWeightOptions {
    bool enable = false;
    float weight_min = 1.0f;
    float weight_max = 5.0f;
    float congestion_alpha = 0.5f;
    int ramp_iterations = 200;
    int update_period = 1;
    float initial_nets_percent = 5.0f;
    float final_nets_percent = 20.0f;
    double slack_norm = 5e-11;
    double length_norm = 100.0;
    float overflow_limit = 0.1f;

    float severity_slack_norm = 2e-10f;
    float severity_ratio_cap = 5.0f;
    float severity_weight_scale = 1.0f;
    float severity_weight_limit = 10.0f;
    float severity_coverage_scale = 1.0f;
    float severity_coverage_limit = 2.0f;

    int top_endpoints = 0;
    float endpoint_slack_threshold = std::numeric_limits<float>::infinity();

    float congestion_gate = 0.9f;
    float hot_bin_fraction = 0.05f;
    float hot_bin_threshold = 0.1f;
  };

  struct CorridorOptions {
    bool enable = false;
    float charge_k = 1.0f;
    int width_bins = 2;
    int top_endpoints = 0;
    int min_path_length = 5;
  };

  struct SpringOptions {
    bool enable = false;
    float k = 1.0f;
    int top_paths = 0;
    float overflow_gate = 0.8f;
    int min_path_length = 5;
  };

  enum class SlackRegime { HighSlack, MediumSlack, LowSlack };

  TimingBase();
  TimingBase(std::shared_ptr<NesterovBaseCommon> nbc,
             rsz::Resizer* rs,
             utl::Logger* log);

  // check whether overflow reached the timingOverflow
  bool isTimingNetWeightOverflow(float overflow);
  void addTimingNetWeightOverflow(int overflow);
  void setTimingNetWeightOverflows(std::vector<int>& overflows);
  void deleteTimingNetWeightOverflow(int overflow);
  void clearTimingNetWeightOverflow();
  size_t getTimingNetWeightOverflowSize() const;

  void setTimingNetWeightMax(float max);
  
  void setDynamicWeightOptions(const DynamicWeightOptions& options);
  void setCorridorOptions(const CorridorOptions& options);
  void setSpringOptions(const SpringOptions& options);
  void setNesterovBases(const std::vector<std::shared_ptr<NesterovBase>>& bases);
  
  void setECPScale(float scale);
  void setECPWeightMax(float max);
  void setECPTopEndpointFrac(float frac);
  
  void updateResizerGuards(const NesterovPlaceVars& vars);

  // updateNetWeight.
  // True: successfully reweighted gnets
  // False: no slacks found
  bool executeTimingDriven(int iter, int total_iter, float avg_overflow, bool run_journal_restore);
  bool executeTimingDriven(bool run_journal_restore);

 private:
  rsz::Resizer* rs_ = nullptr;
  utl::Logger* log_ = nullptr;
  std::shared_ptr<NesterovBaseCommon> nbc_;
  std::vector<std::shared_ptr<NesterovBase>> nbs_;

  std::vector<int> timingNetWeightOverflow_;
  std::vector<int> timingOverflowChk_;
  float net_weight_max_ = 5;
  
  DynamicWeightOptions dynamic_options_;
  CorridorOptions corridor_options_;
  SpringOptions spring_options_;
  
  float ecp_scale_ = 1.0f;
  float ecp_weight_max_ = 2.0f;
  float ecp_top_endpoint_frac_ = 0.05f;

  struct CellLocator { size_t nb_index; size_t cell_index; };
  std::unordered_map<const GCell*, CellLocator> gcell_locator_;
  bool locator_dirty_ = true;

  using PathChain = std::vector<GCell*>;

  void initTimingOverflowChk();
  
  SlackRegime classifyRegime(float wns);
  
  void rebuildGCellLocator();
  std::optional<CellLocator> locateGCell(const GCell* gCell);
  
  std::vector<PathChain> collectPathChains(int top_count, float endpoint_slack_threshold, int min_length);
  
  void applyCorridorMask(const std::vector<PathChain>& chains);
  void clearCorridorMask();
  
  void applyPathSprings(const std::vector<PathChain>& chains);
  void clearPathSprings();
  
  bool crossesCongestedBins(GNet* gNet);
  void applyECPWeights(SlackRegime regime);
  void resetAllNetWeights();
  void applyResizerGuards(SlackRegime regime);

  struct ResizerGuardDefaults {
    double noncrit_budget = 0.0;
    double setup_guard_cap = 0.0;
    double guard_window = 0.0;
    double setup_guard_floor = 0.0;
    bool valid = false;
  };

  ResizerGuardDefaults guard_defaults_;
  bool weights_initialized_ = false;
};

}  // namespace gpl
