// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include "gpl/Replace.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <utility>

#include "AbstractGraphics.h"
#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "graphicsNone.h"
#include "initialPlace.h"
#include "mbff.h"
#include "nesterovBase.h"
#include "nesterovPlace.h"
#include "odb/db.h"
#include "placerBase.h"
#include "routeBase.h"
#include "rsz/Resizer.hh"
#include "sta/StaMain.hh"
#include "timingBase.h"
#include "utl/Logger.h"

namespace gpl {

using utl::GPL;

Replace::Replace(odb::dbDatabase* odb,
                 sta::dbSta* sta,
                 rsz::Resizer* resizer,
                 grt::GlobalRouter* router,
                 utl::Logger* logger)
    : db_(odb), sta_(sta), rs_(resizer), fr_(router), log_(logger)
{
  graphics_ = std::make_unique<GraphicsNone>();
}

Replace::~Replace() = default;

void Replace::setGraphicsInterface(const AbstractGraphics& graphics)
{
  graphics_ = graphics.MakeNew(log_);
}

void Replace::reset()
{
  ip_.reset();
  np_.reset();

  pbc_.reset();
  nbc_.reset();

  pbVec_.clear();
  pbVec_.shrink_to_fit();
  nbVec_.clear();
  nbVec_.shrink_to_fit();

  tb_.reset();
  rb_.reset();

  initialPlaceMaxIter_ = 20;
  initialPlaceMinDiffLength_ = 1500;
  initialPlaceMaxSolverIter_ = 100;
  initialPlaceMaxFanout_ = 200;
  initialPlaceNetWeightScale_ = 800;

  nesterovPlaceMaxIter_ = 5000;
  binGridCntX_ = binGridCntY_ = 0;
  overflow_ = 0.1;
  density_ = 1.0;
  initDensityPenalityFactor_ = 0.00008;
  initWireLengthCoef_ = 0.25;
  minPhiCoef_ = 0.95;
  maxPhiCoef_ = 1.05;
  referenceHpwl_ = 446000000;

  routabilityCheckOverflow_ = 0.3;
  routabilityMaxDensity_ = 0.99;
  routabilityTargetRcMetric_ = 1.01;
  routabilityInflationRatioCoef_ = 2;
  routabilityMaxInflationRatio_ = 3;
  routabilityRcK1_ = routabilityRcK2_ = 1.0;
  routabilityRcK3_ = routabilityRcK4_ = 0.0;
  routabilityMaxInflationIter_ = 4;

  timingDrivenMode_ = true;
  keepResizeBelowOverflow_ = 1.0;
  routabilityDrivenMode_ = true;
  routabilityUseRudy_ = true;
  uniformTargetDensityMode_ = false;
  skipIoMode_ = false;
  disableRevertIfDiverge_ = false;

  padLeft_ = padRight_ = 0;

  timingNetWeightOverflows_.clear();
  timingNetWeightOverflows_.shrink_to_fit();
  timingNetWeightMax_ = 5;

  // Dynamic timing weights
  td_enable_dynamic_weights_ = false;
  td_weight_min_ = 1.0f;
  td_weight_max_ = 5.0f;
  td_congestion_alpha_ = 0.5f;
  td_ramp_iterations_ = 200;
  td_update_period_ = 10;
  td_initial_nets_percent_ = 5.0f;
  td_final_nets_percent_ = 20.0f;
  td_slack_norm_ = 0.1;
  td_length_norm_ = 100.0;
  td_overflow_limit_ = 0.1f;

  // Slack-severity shaping
  td_severity_slack_norm_ = 0.1f;
  td_severity_ratio_cap_ = 5.0f;
  td_severity_weight_scale_ = 1.0f;
  td_severity_weight_limit_ = 10.0f;
  td_severity_coverage_scale_ = 1.0f;
  td_severity_coverage_limit_ = 2.0f;

  // Endpoint-driven selection
  td_top_endpoints_ = 0;
  td_slack_thresh_ = 0.0f;

  // Resizer interaction / guard knobs
  td_noncrit_slack_budget_ns_ = 0.05f;
  td_setup_guard_cap_ns_ = 0.03f;
  td_guard_window_ns_ = 0.0f;
  td_post_cts_hold_floor_ns_ = 0.0f;
  td_rebuffer_clone_gate_fanout_ = 20;
  td_clone_group_fanout_ = 20;
  td_gr_pick_radius_ = 0;

  // Congestion gating & spatial options
  td_congestion_gate_ = 0.9f;
  td_hot_bin_fraction_ = 0.05f;
  td_hot_bin_threshold_ = 0.1f;

  cws_enable_ = false;
  cws_charge_k_ = 1.0f;
  cws_width_bins_ = 2;
  cws_top_endpoints_ = 0;
  cws_min_path_length_ = 5;

  aas_enable_ = false;
  aas_k_ = 1.0f;
  aas_top_paths_ = 0;
  aas_overflow_gate_ = 0.8f;
  aas_min_path_length_ = 5;

  // ECP parameters
  ecp_scale_ = 1.0f;
  ecp_weight_max_ = 2.0f;
  ecp_top_endpoint_frac_ = 0.05f;
}

void Replace::addPlacementCluster(const Cluster& cluster)
{
  clusters_.emplace_back(cluster);
}

void Replace::doIncrementalPlace(int threads)
{
  log_->info(GPL, 6, "Execute incremental mode global placement.");
  if (pbc_ == nullptr) {
    PlacerBaseVars pbVars;
    pbVars.padLeft = padLeft_;
    pbVars.padRight = padRight_;
    pbVars.skipIoMode = skipIoMode_;

    pbc_ = std::make_shared<PlacerBaseCommon>(db_, pbVars, log_);

    pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_));

    for (auto pd : db_->getChip()->getBlock()->getRegions()) {
      for (auto group : pd->getGroups()) {
        pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_, group));
      }
    }

    total_placeable_insts_ = 0;
    for (const auto& pb : pbVec_) {
      total_placeable_insts_ += pb->placeInsts().size();
    }
  }

  // Lock down already placed objects
  int placed_cnt = 0;
  int unplaced_cnt = 0;
  auto block = db_->getChip()->getBlock();
  for (auto inst : block->getInsts()) {
    auto status = inst->getPlacementStatus();
    if (status == odb::dbPlacementStatus::PLACED) {
      pbc_->dbToPb(inst)->lock();
      ++placed_cnt;
    } else if (!status.isPlaced()) {
      ++unplaced_cnt;
    }
  }

  log_->info(GPL, 154, "Identified {} placed instances", placed_cnt);
  log_->info(GPL, 155, "Identified {} not placed instances", unplaced_cnt);

  if (unplaced_cnt == 0) {
    // Everything was already placed so we do the old incremental mode
    // which just skips initial placement and runs nesterov.
    log_->info(GPL,
               156,
               "Identified all instances as placed. Unlocking all instances "
               "and running nesterov from scratch.");
    for (auto& pb : pbVec_) {
      pb->unlockAll();
    }

    doNesterovPlace(threads);
    return;
  }

  // Roughly place the unplaced objects (allow more overflow).
  // Limit iterations to prevent objects drifting too far or
  // non-convergence.
  constexpr float rough_oveflow = 0.2f;
  float previous_overflow = overflow_;
  setTargetOverflow(std::max(rough_oveflow, overflow_));
  doInitialPlace(threads);

  int previous_max_iter = nesterovPlaceMaxIter_;
  initNesterovPlace(threads);
  setNesterovPlaceMaxIter(300);
  int iter = doNesterovPlace(threads);
  setNesterovPlaceMaxIter(previous_max_iter);

  // Finish the overflow resolution from the rough placement
  log_->info(GPL, 133, "Unlocking all instances");
  for (auto& pb : pbVec_) {
    pb->unlockAll();
  }

  setTargetOverflow(previous_overflow);
  if (previous_overflow < rough_oveflow) {
    doNesterovPlace(threads, iter + 1);
  }
}

void Replace::doInitialPlace(int threads)
{
  log_->info(GPL, 5, "Execute conjugate gradient initial placement.");
  if (pbc_ == nullptr) {
    PlacerBaseVars pbVars;
    pbVars.padLeft = padLeft_;
    pbVars.padRight = padRight_;
    pbVars.skipIoMode = skipIoMode_;

    pbc_ = std::make_shared<PlacerBaseCommon>(db_, pbVars, log_);

    pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_));

    for (auto pd : db_->getChip()->getBlock()->getRegions()) {
      for (auto group : pd->getGroups()) {
        pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_, group));
      }
    }

    if (pbVec_.front()->placeInsts().empty()) {
      log_->warn(
          GPL,
          123,
          "No placeable instances in the top-level region. Removing it.");
      pbVec_.erase(pbVec_.begin());
    }

    total_placeable_insts_ = 0;
    for (const auto& pb : pbVec_) {
      total_placeable_insts_ += pb->placeInsts().size();
    }
  }

  InitialPlaceVars ipVars;
  ipVars.maxIter = initialPlaceMaxIter_;
  ipVars.minDiffLength = initialPlaceMinDiffLength_;
  ipVars.maxSolverIter = initialPlaceMaxSolverIter_;
  ipVars.maxFanout = initialPlaceMaxFanout_;
  ipVars.netWeightScale = initialPlaceNetWeightScale_;
  ipVars.debug = gui_debug_initial_;

  std::unique_ptr<InitialPlace> ip(
      new InitialPlace(ipVars, pbc_, pbVec_, graphics_->MakeNew(log_), log_));
  ip_ = std::move(ip);
  ip_->doBicgstabPlace(threads);
}

void Replace::runMBFF(int max_sz,
                      float alpha,
                      float beta,
                      int threads,
                      int num_paths)
{
  MBFF pntset(db_,
              sta_,
              log_,
              rs_,
              threads,
              20,
              num_paths,
              gui_debug_,
              graphics_->MakeNew(log_));
  pntset.Run(max_sz, alpha, beta);
}

bool Replace::initNesterovPlace(int threads)
{
  if (!pbc_) {
    PlacerBaseVars pbVars;
    pbVars.padLeft = padLeft_;
    pbVars.padRight = padRight_;
    pbVars.skipIoMode = skipIoMode_;

    pbc_ = std::make_shared<PlacerBaseCommon>(db_, pbVars, log_);

    pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_));

    for (auto pd : db_->getChip()->getBlock()->getRegions()) {
      for (auto group : pd->getGroups()) {
        pbVec_.push_back(std::make_shared<PlacerBase>(db_, pbc_, log_, group));
      }
    }

    total_placeable_insts_ = 0;
    for (const auto& pb : pbVec_) {
      total_placeable_insts_ += pb->placeInsts().size();
    }
  }

  if (total_placeable_insts_ == 0) {
    log_->warn(GPL, 136, "No placeable instances - skipping placement.");
    return false;
  }

  if (!nbc_) {
    NesterovBaseVars nbVars;
    nbVars.targetDensity = density_;

    if (binGridCntX_ != 0 && binGridCntY_ != 0) {
      nbVars.isSetBinCnt = true;
      nbVars.binCntX = binGridCntX_;
      nbVars.binCntY = binGridCntY_;
      nbVars.minPhiCoef = minPhiCoef_;
      nbVars.maxPhiCoef = maxPhiCoef_;
    }

    nbVars.useUniformTargetDensity = uniformTargetDensityMode_;

    nbc_ = std::make_shared<NesterovBaseCommon>(
        nbVars, pbc_, log_, threads, clusters_);

    for (const auto& pb : pbVec_) {
      nbVec_.push_back(std::make_shared<NesterovBase>(nbVars, pb, nbc_, log_));
    }
  }

  if (!rb_) {
    RouteBaseVars rbVars;
    rbVars.useRudy = routabilityUseRudy_;
    rbVars.maxDensity = routabilityMaxDensity_;
    rbVars.maxInflationIter = routabilityMaxInflationIter_;
    rbVars.targetRC = routabilityTargetRcMetric_;
    rbVars.inflationRatioCoef = routabilityInflationRatioCoef_;
    rbVars.maxInflationRatio = routabilityMaxInflationRatio_;
    rbVars.rcK1 = routabilityRcK1_;
    rbVars.rcK2 = routabilityRcK2_;
    rbVars.rcK3 = routabilityRcK3_;
    rbVars.rcK4 = routabilityRcK4_;

    rb_ = std::make_shared<RouteBase>(rbVars, db_, fr_, nbc_, nbVec_, log_);
  }

  if (!tb_) {
    tb_ = std::make_shared<TimingBase>(nbc_, rs_, log_);
    tb_->setTimingNetWeightOverflows(timingNetWeightOverflows_);
    tb_->setTimingNetWeightMax(timingNetWeightMax_);
  }

  if (!np_) {
    NesterovPlaceVars npVars;

    npVars.referenceHpwl = referenceHpwl_;
    npVars.routability_end_overflow = routabilityCheckOverflow_;
    npVars.keepResizeBelowOverflow = keepResizeBelowOverflow_;
    npVars.initDensityPenalty = initDensityPenalityFactor_;
    npVars.initWireLengthCoef = initWireLengthCoef_;
    npVars.targetOverflow = overflow_;
    npVars.maxNesterovIter = nesterovPlaceMaxIter_;
    npVars.timingDrivenMode = timingDrivenMode_;
    npVars.routability_driven_mode = routabilityDrivenMode_;
    npVars.debug = gui_debug_;
    npVars.debug_pause_iterations = gui_debug_pause_iterations_;
    npVars.debug_update_iterations = gui_debug_update_iterations_;
    npVars.debug_draw_bins = gui_debug_draw_bins_;
    npVars.debug_inst = gui_debug_inst_;
    npVars.debug_start_iter = gui_debug_start_iter_;
    npVars.debug_generate_images = gui_debug_generate_images_;
    npVars.debug_images_path = gui_debug_images_path_;
    npVars.disableRevertIfDiverge = disableRevertIfDiverge_;

    // Dynamic timing weights
    npVars.td_enable_dynamic_weights = td_enable_dynamic_weights_;
    npVars.td_weight_min = td_weight_min_;
    npVars.td_weight_max = td_weight_max_;
    npVars.td_congestion_alpha = td_congestion_alpha_;
    npVars.td_ramp_iterations = td_ramp_iterations_;
    npVars.td_update_period = td_update_period_;
    npVars.td_initial_nets_percent = td_initial_nets_percent_;
    npVars.td_final_nets_percent = td_final_nets_percent_;
    npVars.td_slack_norm = td_slack_norm_;
    npVars.td_length_norm = td_length_norm_;
    npVars.td_overflow_limit = td_overflow_limit_;

    // Slack-severity shaping
    npVars.td_severity_slack_norm = td_severity_slack_norm_;
    npVars.td_severity_ratio_cap = td_severity_ratio_cap_;
    npVars.td_severity_weight_scale = td_severity_weight_scale_;
    npVars.td_severity_weight_limit = td_severity_weight_limit_;
    npVars.td_severity_coverage_scale = td_severity_coverage_scale_;
    npVars.td_severity_coverage_limit = td_severity_coverage_limit_;

    // Endpoint-driven selection
    npVars.td_top_endpoints = td_top_endpoints_;
    npVars.td_slack_thresh = td_slack_thresh_;

    // Resizer interaction / guard knobs
    npVars.td_noncrit_slack_budget_ns = td_noncrit_slack_budget_ns_;
    npVars.td_setup_guard_cap_ns = td_setup_guard_cap_ns_;
    npVars.td_guard_window_ns = td_guard_window_ns_;
    npVars.td_post_cts_hold_floor_ns = td_post_cts_hold_floor_ns_;
    npVars.td_rebuffer_clone_gate_fanout = td_rebuffer_clone_gate_fanout_;
    npVars.td_clone_group_fanout = td_clone_group_fanout_;
    npVars.td_gr_pick_radius = td_gr_pick_radius_;

    // Congestion gating & spatial options
    npVars.td_congestion_gate = td_congestion_gate_;
    npVars.td_hot_bin_fraction = td_hot_bin_fraction_;
    npVars.td_hot_bin_threshold = td_hot_bin_threshold_;

    npVars.cws_enable = cws_enable_;
    npVars.cws_charge_k = cws_charge_k_;
    npVars.cws_width_bins = cws_width_bins_;
    npVars.cws_top_endpoints = cws_top_endpoints_;
    npVars.cws_min_path_length = cws_min_path_length_;

    npVars.aas_enable = aas_enable_;
    npVars.aas_k = aas_k_;
    npVars.aas_top_paths = aas_top_paths_;
    npVars.aas_overflow_gate = aas_overflow_gate_;
    npVars.aas_min_path_length = aas_min_path_length_;

    // ECP parameters
    npVars.ecp_scale = ecp_scale_;
    npVars.ecp_weight_max = ecp_weight_max_;
    npVars.ecp_top_endpoint_frac = ecp_top_endpoint_frac_;

    for (const auto& nb : nbVec_) {
      nb->setNpVars(&npVars);
    }

    np_ = std::make_unique<NesterovPlace>(npVars,
                                          pbc_,
                                          nbc_,
                                          pbVec_,
                                          nbVec_,
                                          rb_,
                                          tb_,
                                          graphics_->MakeNew(log_),
                                          log_);
  }
  return true;
}

int Replace::doNesterovPlace(int threads, int start_iter)
{
  if (!initNesterovPlace(threads)) {
    return 0;
  }

  log_->info(GPL, 7, "Execute nesterov global placement.");
  if (timingDrivenMode_) {
    rs_->resizeSlackPreamble();
  }

  auto start = std::chrono::high_resolution_clock::now();

  int return_do_nesterov = np_->doNesterovPlace(start_iter);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  debugPrint(log_,
             GPL,
             "runtime",
             1,
             "NP->doNesterovPlace() runtime: {} seconds ",
             elapsed.count());

  if (enable_routing_congestion_) {
    fr_->setAllowCongestion(true);
    fr_->setCongestionIterations(0);
    fr_->setCriticalNetsPercentage(0);
    fr_->globalRoute();
  }
  return return_do_nesterov;
}

void Replace::setInitialPlaceMaxIter(int iter)
{
  initialPlaceMaxIter_ = iter;
}

void Replace::setInitialPlaceMinDiffLength(int length)
{
  initialPlaceMinDiffLength_ = length;
}

void Replace::setInitialPlaceMaxSolverIter(int iter)
{
  initialPlaceMaxSolverIter_ = iter;
}

void Replace::setInitialPlaceMaxFanout(int fanout)
{
  initialPlaceMaxFanout_ = fanout;
}

void Replace::setInitialPlaceNetWeightScale(float scale)
{
  initialPlaceNetWeightScale_ = scale;
}

void Replace::setNesterovPlaceMaxIter(int iter)
{
  log_->info(GPL, 158, "Setting nesterov max iterations to {}", iter);
  nesterovPlaceMaxIter_ = iter;
  if (np_) {
    np_->setMaxIters(iter);
  }
}

void Replace::setBinGridCnt(int binGridCntX, int binGridCntY)
{
  binGridCntX_ = binGridCntX;
  binGridCntY_ = binGridCntY;
}

void Replace::setTargetOverflow(float overflow)
{
  log_->info(GPL, 157, "Setting target overflow to {}", overflow);
  overflow_ = overflow;
  if (np_) {
    np_->setTargetOverflow(overflow);
  }
}

void Replace::setTargetDensity(float density)
{
  density_ = density;
}

void Replace::setUniformTargetDensityMode(bool mode)
{
  uniformTargetDensityMode_ = mode;
}

float Replace::getUniformTargetDensity(int threads)
{
  log_->info(GPL, 22, "Initialize gpl and calculate uniform density.");
  log_->redirectStringBegin();

  setSkipIoMode(true);  // in case bterms are not placed

  float density = 1.0f;
  if (initNesterovPlace(threads)) {
    density = nbVec_[0]->getUniformTargetDensity();
  }

  std::string _ = log_->redirectStringEnd();
  return density;
}

void Replace::setInitDensityPenalityFactor(float penaltyFactor)
{
  initDensityPenalityFactor_ = penaltyFactor;
}

void Replace::setInitWireLengthCoef(float coef)
{
  initWireLengthCoef_ = coef;
}

void Replace::setMinPhiCoef(float minPhiCoef)
{
  minPhiCoef_ = minPhiCoef;
}

void Replace::setMaxPhiCoef(float maxPhiCoef)
{
  maxPhiCoef_ = maxPhiCoef;
}

void Replace::setReferenceHpwl(float refHpwl)
{
  referenceHpwl_ = refHpwl;
}

void Replace::setDebug(int pause_iterations,
                       int update_iterations,
                       bool draw_bins,
                       bool initial,
                       odb::dbInst* inst,
                       int start_iter,
                       bool generate_images,
                       std::string images_path)
{
  gui_debug_ = true;
  gui_debug_pause_iterations_ = pause_iterations;
  gui_debug_update_iterations_ = update_iterations;
  gui_debug_draw_bins_ = draw_bins;
  gui_debug_initial_ = initial;
  gui_debug_inst_ = inst;
  gui_debug_start_iter_ = start_iter;
  gui_debug_generate_images_ = generate_images;
  gui_debug_images_path_ = std::move(images_path);
}

void Replace::setDisableRevertIfDiverge(bool mode)
{
  disableRevertIfDiverge_ = mode;
}

void Replace::setEnableRoutingCongestion(bool mode)
{
  enable_routing_congestion_ = mode;
}

void Replace::setSkipIoMode(bool mode)
{
  skipIoMode_ = mode;
}

void Replace::setTimingDrivenMode(bool mode)
{
  timingDrivenMode_ = mode;
}

void Replace::setRoutabilityDrivenMode(bool mode)
{
  routabilityDrivenMode_ = mode;
}

void Replace::setRoutabilityUseGrt(bool mode)
{
  routabilityUseRudy_ = !mode;
}

void Replace::setRoutabilityCheckOverflow(float overflow)
{
  routabilityCheckOverflow_ = overflow;
}

void Replace::setKeepResizeBelowOverflow(float overflow)
{
  keepResizeBelowOverflow_ = overflow;
}

void Replace::setRoutabilityMaxDensity(float density)
{
  routabilityMaxDensity_ = density;
}

void Replace::setRoutabilityMaxInflationIter(int iter)
{
  routabilityMaxInflationIter_ = iter;
}

void Replace::setRoutabilityTargetRcMetric(float rc)
{
  routabilityTargetRcMetric_ = rc;
}

void Replace::setRoutabilityInflationRatioCoef(float coef)
{
  routabilityInflationRatioCoef_ = coef;
}

void Replace::setRoutabilityMaxInflationRatio(float ratio)
{
  routabilityMaxInflationRatio_ = ratio;
}

void Replace::setRoutabilityRcCoefficients(float k1,
                                           float k2,
                                           float k3,
                                           float k4)
{
  routabilityRcK1_ = k1;
  routabilityRcK2_ = k2;
  routabilityRcK3_ = k3;
  routabilityRcK4_ = k4;
}

void Replace::setPadLeft(int pad)
{
  padLeft_ = pad;
}

void Replace::setPadRight(int pad)
{
  padRight_ = pad;
}

void Replace::addTimingNetWeightOverflow(int overflow)
{
  timingNetWeightOverflows_.push_back(overflow);
}

void Replace::setTimingNetWeightMax(float max)
{
  timingNetWeightMax_ = max;
}

void Replace::setTdEnableDynamicWeights(bool enable)
{
  td_enable_dynamic_weights_ = enable;
}

void Replace::setTdWeightMin(float min)
{
  td_weight_min_ = min;
}

void Replace::setTdWeightMax(float max)
{
  td_weight_max_ = max;
}

void Replace::setTdCongestionAlpha(float alpha)
{
  td_congestion_alpha_ = alpha;
}

void Replace::setTdRampIterations(int iterations)
{
  td_ramp_iterations_ = iterations;
}

void Replace::setTdUpdatePeriod(int period)
{
  td_update_period_ = period;
}

void Replace::setTdInitialNetsPercent(float percent)
{
  td_initial_nets_percent_ = percent;
}

void Replace::setTdFinalNetsPercent(float percent)
{
  td_final_nets_percent_ = percent;
}

void Replace::setTdSlackNorm(double norm)
{
  td_slack_norm_ = norm;
}

void Replace::setTdLengthNorm(double norm)
{
  td_length_norm_ = norm;
}

void Replace::setTdOverflowLimit(float limit)
{
  td_overflow_limit_ = limit;
}

void Replace::setTdSeveritySlackNorm(float norm)
{
  td_severity_slack_norm_ = norm;
}

void Replace::setTdSeverityRatioCap(float cap)
{
  td_severity_ratio_cap_ = cap;
}

void Replace::setTdSeverityWeightScale(float scale)
{
  td_severity_weight_scale_ = scale;
}

void Replace::setTdSeverityWeightLimit(float limit)
{
  td_severity_weight_limit_ = limit;
}

void Replace::setTdSeverityCoverageScale(float scale)
{
  td_severity_coverage_scale_ = scale;
}

void Replace::setTdSeverityCoverageLimit(float limit)
{
  td_severity_coverage_limit_ = limit;
}

void Replace::setTdTopEndpoints(int count)
{
  td_top_endpoints_ = count;
}

void Replace::setTdSlackThresh(float thresh)
{
  td_slack_thresh_ = thresh;
}

void Replace::setTdNoncritSlackBudgetNs(float budget)
{
  td_noncrit_slack_budget_ns_ = budget;
}

void Replace::setTdSetupGuardCapNs(float cap)
{
  td_setup_guard_cap_ns_ = cap;
}

void Replace::setTdGuardWindowNs(float window)
{
  td_guard_window_ns_ = window;
}

void Replace::setTdPostCtsHoldFloorNs(float floor)
{
  td_post_cts_hold_floor_ns_ = floor;
}

void Replace::setTdRebufferCloneGateFanout(int fanout)
{
  td_rebuffer_clone_gate_fanout_ = fanout;
}

void Replace::setTdCloneGroupFanout(int fanout)
{
  td_clone_group_fanout_ = fanout;
}

void Replace::setTdGrPickRadius(int radius)
{
  td_gr_pick_radius_ = radius;
}

void Replace::setTdCongestionGate(float gate)
{
  td_congestion_gate_ = gate;
}

void Replace::setTdHotBinFraction(float fraction)
{
  td_hot_bin_fraction_ = fraction;
}

void Replace::setTdHotBinThreshold(float threshold)
{
  td_hot_bin_threshold_ = threshold;
}

void Replace::setCwsEnable(bool enable)
{
  cws_enable_ = enable;
}

void Replace::setCwsChargeK(float k)
{
  cws_charge_k_ = k;
}

void Replace::setCwsWidthBins(int bins)
{
  cws_width_bins_ = bins;
}

void Replace::setCwsTopEndpoints(int count)
{
  cws_top_endpoints_ = count;
}

void Replace::setCwsMinPathLength(int length)
{
  cws_min_path_length_ = length;
}

void Replace::setAasEnable(bool enable)
{
  aas_enable_ = enable;
}

void Replace::setAasK(float k)
{
  aas_k_ = k;
}

void Replace::setAasTopPaths(int count)
{
  aas_top_paths_ = count;
}

void Replace::setAasOverflowGate(float gate)
{
  aas_overflow_gate_ = gate;
}

void Replace::setAasMinPathLength(int length)
{
  aas_min_path_length_ = length;
}

void Replace::setEcpScale(float scale)
{
  ecp_scale_ = scale;
}

void Replace::setEcpWeightMax(float max)
{
  ecp_weight_max_ = max;
}

void Replace::setEcpTopEndpointFrac(float frac)
{
  ecp_top_endpoint_frac_ = frac;
}

}  // namespace gpl
