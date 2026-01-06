// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2020-2025, The OpenROAD Authors

%{
#include "ord/OpenRoad.hh"
#include "gpl/Replace.h"
#include "odb/db.h"

namespace ord {
OpenRoad*
getOpenRoad();

gpl::Replace*
getReplace();

}

using ord::getOpenRoad;
using ord::getReplace;
using gpl::Replace;

%}

%import <std_vector.i>
%import "dbtypes.i"
%import "dbenums.i"
%include "../../Exception.i"

%inline %{

void
placement_cluster_cmd(const std::vector<odb::dbInst*>& cluster)
{
  Replace* replace = getReplace();
  replace->addPlacementCluster(cluster);
}

void 
replace_reset_cmd() 
{
  Replace* replace = getReplace();  
  replace->reset();
}

void 
replace_initial_place_cmd()
{
  Replace* replace = getReplace();
  int threads = ord::OpenRoad::openRoad()->getThreadCount();
  replace->doInitialPlace(threads);
}

void 
replace_nesterov_place_cmd()
{
  Replace* replace = getReplace();
  int threads = ord::OpenRoad::openRoad()->getThreadCount();
  replace->doNesterovPlace(threads);
}


void
replace_run_mbff_cmd(int max_sz, float alpha, float beta, int num_paths) 
{
  Replace* replace = getReplace();
  int threads = ord::OpenRoad::openRoad()->getThreadCount();
  replace->runMBFF(max_sz, alpha, beta, threads, num_paths);   
}


void
set_density_cmd(float density)
{
  Replace* replace = getReplace();
  replace->setTargetDensity(density);
}

void
set_uniform_target_density_mode_cmd(bool uniform)
{
  Replace* replace = getReplace();
  replace->setUniformTargetDensityMode(uniform);
}

void
set_initial_place_max_iter_cmd(int iter)
{
  Replace* replace = getReplace();
  replace->setInitialPlaceMaxIter(iter); 
}

void
set_initial_place_max_fanout_cmd(int fanout)
{
  Replace* replace = getReplace();
  replace->setInitialPlaceMaxFanout(fanout);
}

void
set_nesv_place_iter_cmd(int iter)
{
  Replace* replace = getReplace();
  replace->setNesterovPlaceMaxIter(iter);
}

void
set_bin_grid_cnt_cmd(int cnt_x, int cnt_y)
{
  Replace* replace = getReplace();
  replace->setBinGridCnt(cnt_x, cnt_y);
}

void
set_overflow_cmd(float overflow)
{
  Replace* replace = getReplace();
  replace->setTargetOverflow(overflow);
}

void
set_min_phi_coef_cmd(float min_phi_coef)
{
  Replace* replace = getReplace();
  replace->setMinPhiCoef(min_phi_coef);
}

void
set_max_phi_coef_cmd(float max_phi_coef) 
{
  Replace* replace = getReplace();
  replace->setMaxPhiCoef(max_phi_coef);
}

void
set_reference_hpwl_cmd(float reference_hpwl)
{
  Replace* replace = getReplace();
  replace->setReferenceHpwl(reference_hpwl);
}

void
set_init_density_penalty_factor_cmd(float penaltyFactor)
{
  Replace* replace = getReplace();
  replace->setInitDensityPenalityFactor(penaltyFactor);
}

void
set_init_wirelength_coef_cmd(float coef)
{
  Replace* replace = getReplace();
  replace->setInitWireLengthCoef(coef);
}

void
replace_incremental_place_cmd()
{
  Replace* replace = getReplace();
  int threads = ord::OpenRoad::openRoad()->getThreadCount();
  replace->doIncrementalPlace(threads);
}


void set_timing_driven_mode(bool timing_driven)
{
  Replace* replace = getReplace();
  replace->setTimingDrivenMode(timing_driven);
}


void
set_keep_resize_below_overflow_cmd(float overflow) 
{
  Replace* replace = getReplace();
  replace->setKeepResizeBelowOverflow(overflow);
}

void
set_routability_driven_mode(bool routability_driven)
{
  Replace* replace = getReplace();
  replace->setRoutabilityDrivenMode(routability_driven);
}

void
set_routability_use_grt(bool use_grt)
{
  Replace* replace = getReplace();
  replace->setRoutabilityUseGrt(use_grt);
}

void
set_routability_check_overflow_cmd(float overflow) 
{
  Replace* replace = getReplace();
  replace->setRoutabilityCheckOverflow(overflow);
}
 
void
set_routability_max_density_cmd(float density) 
{
  Replace* replace = getReplace();
  replace->setRoutabilityMaxDensity(density);
}

void
set_routability_max_inflation_iter_cmd(int iter) 
{
  Replace* replace = getReplace();
  replace->setRoutabilityMaxInflationIter(iter);
}

void
set_routability_target_rc_metric_cmd(float rc)
{
  Replace* replace = getReplace();
  replace->setRoutabilityTargetRcMetric(rc);
}

void
set_routability_inflation_ratio_coef_cmd(float coef)
{
  Replace* replace = getReplace();
  replace->setRoutabilityInflationRatioCoef(coef);
}

void
set_routability_max_inflation_ratio_cmd(float ratio) 
{
  Replace* replace = getReplace();
  replace->setRoutabilityMaxInflationRatio(ratio);
}

void
set_routability_rc_coefficients_cmd(float k1,
                                    float k2,
                                    float k3,
                                    float k4)
{
  Replace* replace = getReplace();
  replace->setRoutabilityRcCoefficients(k1, k2, k3, k4);
}


void
set_pad_left_cmd(int pad) 
{
  Replace* replace = getReplace();
  replace->setPadLeft(pad);
}

void
set_pad_right_cmd(int pad) 
{
  Replace* replace = getReplace();
  replace->setPadRight(pad);
}

void
set_skip_io_mode_cmd(bool mode) 
{
  Replace* replace = getReplace();
  replace->setSkipIoMode(mode);
}

void
set_disable_revert_if_diverge(bool disable_revert_if_diverge)
{
  Replace* replace = getReplace();
  replace->setDisableRevertIfDiverge(disable_revert_if_diverge);
}

void
set_enable_routing_congestion(bool enable_routing_congestion)
{
  Replace* replace = getReplace();
  replace->setEnableRoutingCongestion(enable_routing_congestion);
}

float
get_global_placement_uniform_density_cmd() 
{
  Replace* replace = getReplace();
  int threads = ord::OpenRoad::openRoad()->getThreadCount();
  return replace->getUniformTargetDensity(threads);
}

void 
add_timing_net_reweight_overflow_cmd(int overflow)
{
  Replace* replace = getReplace();
  return replace->addTimingNetWeightOverflow(overflow);
}

void
set_timing_driven_net_weight_max_cmd(float max)
{
  Replace* replace = getReplace();
  return replace->setTimingNetWeightMax(max);
}

void set_td_enable_dynamic_weights_cmd(bool enable) { getReplace()->setTdEnableDynamicWeights(enable); }
void set_td_weight_min_cmd(float min) { getReplace()->setTdWeightMin(min); }
void set_td_weight_max_cmd(float max) { getReplace()->setTdWeightMax(max); }
void set_td_congestion_alpha_cmd(float alpha) { getReplace()->setTdCongestionAlpha(alpha); }
void set_td_ramp_iterations_cmd(int iter) { getReplace()->setTdRampIterations(iter); }
void set_td_update_period_cmd(int period) { getReplace()->setTdUpdatePeriod(period); }
void set_td_initial_nets_percent_cmd(float percent) { getReplace()->setTdInitialNetsPercent(percent); }
void set_td_final_nets_percent_cmd(float percent) { getReplace()->setTdFinalNetsPercent(percent); }
void set_td_slack_norm_cmd(double norm) { getReplace()->setTdSlackNorm(norm); }
void set_td_length_norm_cmd(double norm) { getReplace()->setTdLengthNorm(norm); }
void set_td_overflow_limit_cmd(float limit) { getReplace()->setTdOverflowLimit(limit); }
void set_td_severity_slack_norm_cmd(float norm) { getReplace()->setTdSeveritySlackNorm(norm); }
void set_td_severity_ratio_cap_cmd(float cap) { getReplace()->setTdSeverityRatioCap(cap); }
void set_td_severity_weight_scale_cmd(float scale) { getReplace()->setTdSeverityWeightScale(scale); }
void set_td_severity_weight_limit_cmd(float limit) { getReplace()->setTdSeverityWeightLimit(limit); }
void set_td_severity_coverage_scale_cmd(float scale) { getReplace()->setTdSeverityCoverageScale(scale); }
void set_td_severity_coverage_limit_cmd(float limit) { getReplace()->setTdSeverityCoverageLimit(limit); }
void set_td_top_endpoints_cmd(int count) { getReplace()->setTdTopEndpoints(count); }
void set_td_slack_thresh_cmd(float thresh) { getReplace()->setTdSlackThresh(thresh); }
void set_td_noncrit_slack_budget_ns_cmd(float budget) { getReplace()->setTdNoncritSlackBudgetNs(budget); }
void set_td_setup_guard_cap_ns_cmd(float cap) { getReplace()->setTdSetupGuardCapNs(cap); }
void set_td_guard_window_ns_cmd(float window) { getReplace()->setTdGuardWindowNs(window); }
void set_td_post_cts_hold_floor_ns_cmd(float floor) { getReplace()->setTdPostCtsHoldFloorNs(floor); }
void set_td_rebuffer_clone_gate_fanout_cmd(int fanout) { getReplace()->setTdRebufferCloneGateFanout(fanout); }
void set_td_clone_group_fanout_cmd(int fanout) { getReplace()->setTdCloneGroupFanout(fanout); }
void set_td_gr_pick_radius_cmd(int radius) { getReplace()->setTdGrPickRadius(radius); }
void set_td_congestion_gate_cmd(float gate) { getReplace()->setTdCongestionGate(gate); }
void set_td_hot_bin_fraction_cmd(float frac) { getReplace()->setTdHotBinFraction(frac); }
void set_td_hot_bin_threshold_cmd(float thresh) { getReplace()->setTdHotBinThreshold(thresh); }
void set_cws_enable_cmd(bool enable) { getReplace()->setCwsEnable(enable); }
void set_cws_charge_k_cmd(float k) { getReplace()->setCwsChargeK(k); }
void set_cws_width_bins_cmd(int bins) { getReplace()->setCwsWidthBins(bins); }
void set_cws_top_endpoints_cmd(int count) { getReplace()->setCwsTopEndpoints(count); }
void set_cws_min_path_length_cmd(int len) { getReplace()->setCwsMinPathLength(len); }
void set_aas_enable_cmd(bool enable) { getReplace()->setAasEnable(enable); }
void set_aas_k_cmd(float k) { getReplace()->setAasK(k); }
void set_aas_top_paths_cmd(int count) { getReplace()->setAasTopPaths(count); }
void set_aas_overflow_gate_cmd(float gate) { getReplace()->setAasOverflowGate(gate); }
void set_aas_min_path_length_cmd(int len) { getReplace()->setAasMinPathLength(len); }
void set_ecp_scale_cmd(float scale) { getReplace()->setEcpScale(scale); }
void set_ecp_weight_max_cmd(float max) { getReplace()->setEcpWeightMax(max); }
void set_ecp_top_endpoint_frac_cmd(float frac) { getReplace()->setEcpTopEndpointFrac(frac); }

void
set_debug_cmd(int pause_iterations,
              int update_iterations,
              bool draw_bins,
              bool initial,
              const char* inst_name,
              int start_iter,
              bool generate_images,
              const char* images_path)
{
  Replace* replace = getReplace();
  odb::dbInst* inst = nullptr;
  if (inst_name) {
    auto block = ord::OpenRoad::openRoad()->getDb()->getChip()->getBlock();
    inst = block->findInst(inst_name);
  }

  std::string resolved_path = (images_path && *images_path)
                                  ? images_path
                                  : "REPORTS_DIR";

  replace->setDebug(pause_iterations, update_iterations, draw_bins,
                    initial, inst, start_iter, generate_images,
                    resolved_path);
}

%} // inline
