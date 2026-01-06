# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2018-2025, The OpenROAD Authors

sta::define_cmd_args "global_placement" {\
    [-skip_initial_place]\
    [-skip_nesterov_place]\
    [-timing_driven]\
    [-routability_driven]\
    [-disable_timing_driven]\
    [-disable_routability_driven]\
    [-incremental]\
    [-skip_io]\
    [-bin_grid_count grid_count]\
    [-density target_density]\
    [-init_density_penalty init_density_penalty]\
    [-init_wirelength_coef init_wirelength_coef]\
    [-min_phi_coef min_phi_coef]\
    [-max_phi_coef max_phi_coef]\
    [-reference_hpwl reference_hpwl]\
    [-overflow overflow]\
    [-initial_place_max_iter initial_place_max_iter]\
    [-initial_place_max_fanout initial_place_max_fanout]\
    [-routability_use_grt]\
    [-routability_target_rc_metric routability_target_rc_metric]\
    [-routability_check_overflow routability_check_overflow]\
    [-routability_max_density routability_max_density]\
    [-routability_max_inflation_iter routability_max_inflation_iter]\
    [-routability_inflation_ratio_coef routability_inflation_ratio_coef]\
    [-routability_max_inflation_ratio routability_max_inflation_ratio]\
    [-routability_rc_coefficients routability_rc_coefficients]\
    [-keep_resize_below_overflow keep_resize_below_overflow]\
    [-timing_driven_net_reweight_overflow timing_driven_net_reweight_overflow]\
    [-timing_driven_net_weight_max timing_driven_net_weight_max]\
    [-timing_driven_nets_percentage timing_driven_nets_percentage]\
    [-td_enable_dynamic_weights]\
    [-td_weight_min td_weight_min]\
    [-td_weight_max td_weight_max]\
    [-td_congestion_alpha td_congestion_alpha]\
    [-td_ramp_iterations td_ramp_iterations]\
    [-td_update_period td_update_period]\
    [-td_initial_nets_percent td_initial_nets_percent]\
    [-td_final_nets_percent td_final_nets_percent]\
    [-td_slack_norm td_slack_norm]\
    [-td_length_norm td_length_norm]\
    [-td_overflow_limit td_overflow_limit]\
    [-td_severity_slack_norm td_severity_slack_norm]\
    [-td_severity_ratio_cap td_severity_ratio_cap]\
    [-td_severity_weight_scale td_severity_weight_scale]\
    [-td_severity_weight_limit td_severity_weight_limit]\
    [-td_severity_coverage_scale td_severity_coverage_scale]\
    [-td_severity_coverage_limit td_severity_coverage_limit]\
    [-td_top_endpoints td_top_endpoints]\
    [-td_slack_thresh td_slack_thresh]\
    [-td_noncrit_slack_budget_ns td_noncrit_slack_budget_ns]\
    [-td_setup_guard_cap_ns td_setup_guard_cap_ns]\
    [-td_guard_window_ns td_guard_window_ns]\
    [-td_post_cts_hold_floor_ns td_post_cts_hold_floor_ns]\
    [-td_rebuffer_clone_gate_fanout td_rebuffer_clone_gate_fanout]\
    [-td_clone_group_fanout td_clone_group_fanout]\
    [-td_gr_pick_radius td_gr_pick_radius]\
    [-td_congestion_gate td_congestion_gate]\
    [-td_hot_bin_fraction td_hot_bin_fraction]\
    [-td_hot_bin_threshold td_hot_bin_threshold]\
    [-cws_enable]\
    [-cws_charge_k cws_charge_k]\
    [-cws_width_bins cws_width_bins]\
    [-cws_top_endpoints cws_top_endpoints]\
    [-cws_min_path_length cws_min_path_length]\
    [-aas_enable]\
    [-aas_k aas_k]\
    [-aas_top_paths aas_top_paths]\
    [-aas_overflow_gate aas_overflow_gate]\
    [-aas_min_path_length aas_min_path_length]\
    [-ecp_scale ecp_scale]\
    [-ecp_weight_max ecp_weight_max]\
    [-ecp_top_endpoint_frac ecp_top_endpoint_frac]\
    [-pad_left pad_left]\
    [-pad_right pad_right]\
    [-disable_revert_if_diverge]\
    [-enable_routing_congestion]
}

proc global_placement { args } {
  sta::parse_key_args "global_placement" args \
    keys {-bin_grid_count -density \
      -init_density_penalty -init_wirelength_coef \
      -min_phi_coef -max_phi_coef -overflow \
      -reference_hpwl \
      -initial_place_max_iter -initial_place_max_fanout \
      -routability_check_overflow -routability_max_density \
      -routability_max_inflation_iter \
      -routability_target_rc_metric \
      -routability_inflation_ratio_coef \
      -routability_max_inflation_ratio \
      -routability_rc_coefficients \
      -timing_driven_net_reweight_overflow \
      -timing_driven_net_weight_max \
      -timing_driven_nets_percentage \
      -keep_resize_below_overflow \
      -td_weight_min -td_weight_max -td_congestion_alpha \
      -td_ramp_iterations -td_update_period \
      -td_initial_nets_percent -td_final_nets_percent \
      -td_slack_norm -td_length_norm -td_overflow_limit \
      -td_severity_slack_norm -td_severity_ratio_cap \
      -td_severity_weight_scale -td_severity_weight_limit \
      -td_severity_coverage_scale -td_severity_coverage_limit \
      -td_top_endpoints -td_slack_thresh \
      -td_noncrit_slack_budget_ns -td_setup_guard_cap_ns \
      -td_guard_window_ns -td_post_cts_hold_floor_ns \
      -td_rebuffer_clone_gate_fanout -td_clone_group_fanout \
      -td_gr_pick_radius -td_congestion_gate \
      -td_hot_bin_fraction -td_hot_bin_threshold \
      -cws_charge_k -cws_width_bins -cws_top_endpoints -cws_min_path_length \
      -aas_k -aas_top_paths -aas_overflow_gate -aas_min_path_length \
      -ecp_scale -ecp_weight_max -ecp_top_endpoint_frac \
      -pad_left -pad_right} \
    flags {-skip_initial_place \
      -skip_nesterov_place \
      -timing_driven \
      -routability_driven \
      -routability_use_grt \
      -disable_timing_driven \
      -disable_routability_driven \
      -skip_io \
      -incremental \
      -disable_revert_if_diverge \
      -enable_routing_congestion \
      -td_enable_dynamic_weights \
      -cws_enable -aas_enable}

  # flow control for the initial placement phase
  if { [info exists flags(-skip_initial_place)] } {
    gpl::set_initial_place_max_iter_cmd 0
  } elseif { [info exists keys(-initial_place_max_iter)] } {
    set initial_place_max_iter $keys(-initial_place_max_iter)
    sta::check_positive_integer "-initial_place_max_iter" $initial_place_max_iter
    gpl::set_initial_place_max_iter_cmd $initial_place_max_iter
  }

  # Skip IO handling
  set skip_io [info exists flags(-skip_io)]
  gpl::set_skip_io_mode_cmd $skip_io
  if { $skip_io } {
    # no IO legalization run during the GP skip-io stage
    gpl::set_initial_place_max_iter_cmd 0
  }

  set timing_driven [info exists flags(-timing_driven)]
  gpl::set_timing_driven_mode $timing_driven
  if { $timing_driven } {
    if { $skip_io } {
      utl::warn "GPL" 150 "-skip_io will disable timing driven mode."
      gpl::set_timing_driven_mode 0
    }

    if { [info exists flags(-td_enable_dynamic_weights)] } { gpl::set_td_enable_dynamic_weights_cmd 1 }
    if { [info exists keys(-td_weight_min)] } { gpl::set_td_weight_min_cmd $keys(-td_weight_min) }
    if { [info exists keys(-td_weight_max)] } { gpl::set_td_weight_max_cmd $keys(-td_weight_max) }
    if { [info exists keys(-td_congestion_alpha)] } { gpl::set_td_congestion_alpha_cmd $keys(-td_congestion_alpha) }
    if { [info exists keys(-td_ramp_iterations)] } { gpl::set_td_ramp_iterations_cmd $keys(-td_ramp_iterations) }
    if { [info exists keys(-td_update_period)] } { gpl::set_td_update_period_cmd $keys(-td_update_period) }
    if { [info exists keys(-td_initial_nets_percent)] } { gpl::set_td_initial_nets_percent_cmd $keys(-td_initial_nets_percent) }
    if { [info exists keys(-td_final_nets_percent)] } { gpl::set_td_final_nets_percent_cmd $keys(-td_final_nets_percent) }
    if { [info exists keys(-td_slack_norm)] } { gpl::set_td_slack_norm_cmd $keys(-td_slack_norm) }
    if { [info exists keys(-td_length_norm)] } { gpl::set_td_length_norm_cmd $keys(-td_length_norm) }
    if { [info exists keys(-td_overflow_limit)] } { gpl::set_td_overflow_limit_cmd $keys(-td_overflow_limit) }
    
    if { [info exists keys(-td_severity_slack_norm)] } { gpl::set_td_severity_slack_norm_cmd $keys(-td_severity_slack_norm) }
    if { [info exists keys(-td_severity_ratio_cap)] } { gpl::set_td_severity_ratio_cap_cmd $keys(-td_severity_ratio_cap) }
    if { [info exists keys(-td_severity_weight_scale)] } { gpl::set_td_severity_weight_scale_cmd $keys(-td_severity_weight_scale) }
    if { [info exists keys(-td_severity_weight_limit)] } { gpl::set_td_severity_weight_limit_cmd $keys(-td_severity_weight_limit) }
    if { [info exists keys(-td_severity_coverage_scale)] } { gpl::set_td_severity_coverage_scale_cmd $keys(-td_severity_coverage_scale) }
    if { [info exists keys(-td_severity_coverage_limit)] } { gpl::set_td_severity_coverage_limit_cmd $keys(-td_severity_coverage_limit) }
    
    if { [info exists keys(-td_top_endpoints)] } { gpl::set_td_top_endpoints_cmd $keys(-td_top_endpoints) }
    if { [info exists keys(-td_slack_thresh)] } { gpl::set_td_slack_thresh_cmd $keys(-td_slack_thresh) }
    
    if { [info exists keys(-td_noncrit_slack_budget_ns)] } { gpl::set_td_noncrit_slack_budget_ns_cmd $keys(-td_noncrit_slack_budget_ns) }
    if { [info exists keys(-td_setup_guard_cap_ns)] } { gpl::set_td_setup_guard_cap_ns_cmd $keys(-td_setup_guard_cap_ns) }
    if { [info exists keys(-td_guard_window_ns)] } { gpl::set_td_guard_window_ns_cmd $keys(-td_guard_window_ns) }
    if { [info exists keys(-td_post_cts_hold_floor_ns)] } { gpl::set_td_post_cts_hold_floor_ns_cmd $keys(-td_post_cts_hold_floor_ns) }
    if { [info exists keys(-td_rebuffer_clone_gate_fanout)] } { gpl::set_td_rebuffer_clone_gate_fanout_cmd $keys(-td_rebuffer_clone_gate_fanout) }
    if { [info exists keys(-td_clone_group_fanout)] } { gpl::set_td_clone_group_fanout_cmd $keys(-td_clone_group_fanout) }
    if { [info exists keys(-td_gr_pick_radius)] } { gpl::set_td_gr_pick_radius_cmd $keys(-td_gr_pick_radius) }
    
    if { [info exists keys(-td_congestion_gate)] } { gpl::set_td_congestion_gate_cmd $keys(-td_congestion_gate) }
    if { [info exists keys(-td_hot_bin_fraction)] } { gpl::set_td_hot_bin_fraction_cmd $keys(-td_hot_bin_fraction) }
    if { [info exists keys(-td_hot_bin_threshold)] } { gpl::set_td_hot_bin_threshold_cmd $keys(-td_hot_bin_threshold) }
    
    if { [info exists flags(-cws_enable)] } { gpl::set_cws_enable_cmd 1 }
    if { [info exists keys(-cws_charge_k)] } { gpl::set_cws_charge_k_cmd $keys(-cws_charge_k) }
    if { [info exists keys(-cws_width_bins)] } { gpl::set_cws_width_bins_cmd $keys(-cws_width_bins) }
    if { [info exists keys(-cws_top_endpoints)] } { gpl::set_cws_top_endpoints_cmd $keys(-cws_top_endpoints) }
    if { [info exists keys(-cws_min_path_length)] } { gpl::set_cws_min_path_length_cmd $keys(-cws_min_path_length) }
    
    if { [info exists flags(-aas_enable)] } { gpl::set_aas_enable_cmd 1 }
    if { [info exists keys(-aas_k)] } { gpl::set_aas_k_cmd $keys(-aas_k) }
    if { [info exists keys(-aas_top_paths)] } { gpl::set_aas_top_paths_cmd $keys(-aas_top_paths) }
    if { [info exists keys(-aas_overflow_gate)] } { gpl::set_aas_overflow_gate_cmd $keys(-aas_overflow_gate) }
    if { [info exists keys(-aas_min_path_length)] } { gpl::set_aas_min_path_length_cmd $keys(-aas_min_path_length) }
    
    if { [info exists keys(-ecp_scale)] } { gpl::set_ecp_scale_cmd $keys(-ecp_scale) }
    if { [info exists keys(-ecp_weight_max)] } { gpl::set_ecp_weight_max_cmd $keys(-ecp_weight_max) }
    if { [info exists keys(-ecp_top_endpoint_frac)] } { gpl::set_ecp_top_endpoint_frac_cmd $keys(-ecp_top_endpoint_frac) }

    if { [info exists keys(-timing_driven_net_reweight_overflow)] } {
      set overflow_list $keys(-timing_driven_net_reweight_overflow)
    } else {
      set overflow_list [list 64 20]
    }

    foreach overflow $overflow_list {
      gpl::add_timing_net_reweight_overflow_cmd $overflow
    }

    # timing driven check overflow to keep resizer changes (non-virtual resizer)
    if { [info exists keys(-keep_resize_below_overflow)] } {
      set keep_resize_below_overflow $keys(-keep_resize_below_overflow)
      sta::check_positive_float "-keep_resize_below_overflow" $keep_resize_below_overflow
      gpl::set_keep_resize_below_overflow_cmd $keep_resize_below_overflow
    }

    if { [info exists keys(-timing_driven_net_weight_max)] } {
      gpl::set_timing_driven_net_weight_max_cmd $keys(-timing_driven_net_weight_max)
    }

    if { [info exists keys(-timing_driven_nets_percentage)] } {
      rsz::set_worst_slack_nets_percent $keys(-timing_driven_nets_percentage)
    }
  }

  if { [info exists flags(-disable_timing_driven)] } {
    utl::warn "GPL" 115 "-disable_timing_driven is deprecated."
  }

  set routability_driven [info exists flags(-routability_driven)]
  gpl::set_routability_driven_mode $routability_driven
  if { $routability_driven } {
    if { $skip_io } {
      utl::warn "GPL" 151 "-skip_io will disable routability driven mode."
      gpl::set_routability_driven_mode 0
    }
  }
  if { [info exists flags(-disable_routability_driven)] } {
    utl::warn "GPL" 116 "-disable_routability_driven is deprecated."
  }

  set routability_use_grt [info exists flags(-routability_use_grt)]
  gpl::set_routability_use_grt $routability_use_grt
  if { $routability_driven } {
    if { $routability_use_grt } {
      utl::warn "GPL" 152 \
        "Using GRT FastRoute instead of default RUDY for congestion in routability driven."
    }
  }

  # Disable revert to saved snapshot if a divergence is detected.
  set disable_revert_if_diverge [info exists flags(-disable_revert_if_diverge)]
  gpl::set_disable_revert_if_diverge $disable_revert_if_diverge
  if { $disable_revert_if_diverge } {
    utl::info "GPL" 153 \
      "Revert-to-snapshot on divergence detection is disabled."
  }

  set enable_routing_congestion [info exists flags(-enable_routing_congestion)]
  gpl::set_enable_routing_congestion $enable_routing_congestion

  if { [info exists keys(-initial_place_max_fanout)] } {
    set initial_place_max_fanout $keys(-initial_place_max_fanout)
    sta::check_positive_integer "-initial_place_max_fanout" $initial_place_max_fanout
    gpl::set_initial_place_max_fanout_cmd $initial_place_max_fanout
  }

  # density settings
  set target_density 0.7
  set uniform_mode 0

  if { [info exists keys(-density)] } {
    set target_density $keys(-density)
  }

  if { $target_density == "uniform" } {
    set uniform_mode 1
  } else {
    sta::check_positive_float "-density" $target_density
    if { $target_density > 1.0 } {
      utl::error GPL 135 "Target density must be in \[0, 1\]."
    }
    gpl::set_density_cmd $target_density
  }

  gpl::set_uniform_target_density_mode_cmd $uniform_mode

  if { [info exists keys(-routability_max_density)] } {
    set routability_max_density $keys(-routability_max_density)
    sta::check_positive_float "-routability_max_density" $routability_max_density
    gpl::set_routability_max_density_cmd $routability_max_density
  }


  # parameter to control the RePlAce divergence
  if { [info exists keys(-min_phi_coef)] } {
    set min_phi_coef $keys(-min_phi_coef)
    sta::check_positive_float "-min_phi_coef" $min_phi_coef
    gpl::set_min_phi_coef_cmd $min_phi_coef
  }

  if { [info exists keys(-max_phi_coef)] } {
    set max_phi_coef $keys(-max_phi_coef)
    sta::check_positive_float "-max_phi_coef" $max_phi_coef
    gpl::set_max_phi_coef_cmd $max_phi_coef
  }

  if { [info exists keys(-init_density_penalty)] } {
    set density_penalty $keys(-init_density_penalty)
    sta::check_positive_float "-init_density_penalty" $density_penalty
    gpl::set_init_density_penalty_factor_cmd $density_penalty
  }

  if { [info exists keys(-init_wirelength_coef)] } {
    set coef $keys(-init_wirelength_coef)
    sta::check_positive_float "-init_wirelength_coef" $coef
    gpl::set_init_wirelength_coef_cmd $coef
  }

  if { [info exists keys(-reference_hpwl)] } {
    set reference_hpwl $keys(-reference_hpwl)
    sta::check_positive_float "-reference_hpwl" $reference_hpwl
    gpl::set_reference_hpwl_cmd $reference_hpwl
  }


  if { [info exists keys(-bin_grid_count)] } {
    set bin_grid_count $keys(-bin_grid_count)
    sta::check_positive_integer "-bin_grid_count" $bin_grid_count
    gpl::set_bin_grid_cnt_cmd $bin_grid_count $bin_grid_count
  }

  # overflow
  if { [info exists keys(-overflow)] } {
    set overflow $keys(-overflow)
    sta::check_positive_float "-overflow" $overflow
    gpl::set_overflow_cmd $overflow
  }

  # routability check overflow
  if { [info exists keys(-routability_check_overflow)] } {
    set routability_check_overflow $keys(-routability_check_overflow)
    sta::check_positive_float "-routability_check_overflow" $routability_check_overflow
    gpl::set_routability_check_overflow_cmd $routability_check_overflow
  }

  # routability inflation iter
  if { [info exists keys(-routability_max_inflation_iter)] } {
    set routability_max_inflation_iter $keys(-routability_max_inflation_iter)
    sta::check_positive_float "-routability_max_inflation_iter" $routability_max_inflation_iter
    gpl::set_routability_max_inflation_iter_cmd $routability_max_inflation_iter
  }

  # routability inflation iter
  if { [info exists keys(-routability_target_rc_metric)] } {
    set target_rc_metric $keys(-routability_target_rc_metric)
    sta::check_positive_float "-routability_target_rc_metric" $target_rc_metric
    gpl::set_routability_target_rc_metric_cmd $target_rc_metric
  }

  # routability inflation ratio coef
  if { [info exists keys(-routability_inflation_ratio_coef)] } {
    set ratio_coef $keys(-routability_inflation_ratio_coef)
    sta::check_positive_float "-routability_inflation_ratio_coef" $ratio_coef
    gpl::set_routability_inflation_ratio_coef_cmd $ratio_coef
  }

  # routability max inflation ratio
  if { [info exists keys(-routability_max_inflation_ratio)] } {
    set max_inflation_ratio $keys(-routability_max_inflation_ratio)
    sta::check_positive_float "-routability_max_inflation_ratio" $max_inflation_ratio
    gpl::set_routability_max_inflation_ratio_cmd $max_inflation_ratio
  }

  # routability rc coefficients control
  if { [info exists keys(-routability_rc_coefficients)] } {
    set rc_coefficients $keys(-routability_rc_coefficients)
    set k1 [lindex $rc_coefficients 0]
    set k2 [lindex $rc_coefficients 1]
    set k3 [lindex $rc_coefficients 2]
    set k4 [lindex $rc_coefficients 3]
    gpl::set_routability_rc_coefficients_cmd $k1 $k2 $k3 $k4
  }

  # temp code.
  if { [info exists keys(-pad_left)] } {
    set pad_left $keys(-pad_left)
    sta::check_positive_integer "-pad_left" $pad_left
    gpl::set_pad_left_cmd $pad_left
  }
  if { [info exists keys(-pad_right)] } {
    set pad_right $keys(-pad_right)
    sta::check_positive_integer "-pad_right" $pad_right
    gpl::set_pad_right_cmd $pad_right
  }

  if { [ord::db_has_rows] } {
    sta::check_argc_eq0 "global_placement" $args

    if { [info exists flags(-incremental)] } {
      gpl::replace_incremental_place_cmd
    } else {
      gpl::replace_initial_place_cmd

      if { ![info exists flags(-skip_nesterov_place)] } {
        gpl::replace_nesterov_place_cmd
      }
    }
    gpl::replace_reset_cmd
  } else {
    utl::error GPL 130 "No rows defined in design. Use initialize_floorplan to add rows."
  }
}


sta::define_cmd_args "cluster_flops" {\
    [-tray_weight tray_weight]\
    [-timing_weight timing_weight]\
    [-max_split_size max_split_size]\
    [-num_paths num_paths]\
}

proc cluster_flops { args } {
  sta::parse_key_args "cluster_flops" args \
    keys { -tray_weight -timing_weight -max_split_size -num_paths } \
    flags {}

  if { [ord::get_db_block] == "NULL" } {
    utl::error GPL 113 "No design block found."
  }

  set tray_weight 32.0
  set timing_weight 0.1
  set max_split_size 500
  set num_paths 0

  if { [info exists keys(-tray_weight)] } {
    set tray_weight $keys(-tray_weight)
  }

  if { [info exists keys(-timing_weight)] } {
    set timing_weight $keys(-timing_weight)
  }

  if { [info exists keys(-max_split_size)] } {
    set max_split_size $keys(-max_split_size)
  }

  if { [info exists keys(-num_paths)] } {
    set num_paths $keys(-num_paths)
  }

  gpl::replace_run_mbff_cmd $max_split_size $tray_weight $timing_weight $num_paths
}

proc global_placement_debug { args } {
  sta::parse_key_args "global_placement_debug" args \
    keys {-pause -update -inst -start_iter -images_path} \
    flags {-draw_bins -initial -generate_images} ;# checker off

  if { [ord::get_db_block] == "NULL" } {
    utl::error GPL 117 "No design block found."
  }

  set pause 10
  if { [info exists keys(-pause)] } {
    set pause $keys(-pause)
    sta::check_positive_integer "-pause" $pause
  }

  set update 10
  if { [info exists keys(-update)] } {
    set update $keys(-update)
    sta::check_positive_integer "-update" $update
  }

  set inst ""
  if { [info exists keys(-inst)] } {
    set inst $keys(-inst)
  }

  set start_iter 0
  if { [info exists keys(-start_iter)] } {
    set start_iter $keys(-start_iter)
    sta::check_positive_integer "-start_iter" $start_iter
  }

  set draw_bins [info exists flags(-draw_bins)]
  set initial [info exists flags(-initial)]
  set generate_images [info exists flags(-generate_images)]

  set images_path ""
  if { [info exists keys(-images_path)] } {
    set images_path $keys(-images_path)
  }

  gpl::set_debug_cmd $pause $update $draw_bins $initial \
    $inst $start_iter $generate_images $images_path
}

sta::define_cmd_args "placement_cluster" {}

proc placement_cluster { args } {
  sta::parse_key_args "placement_cluster" args \
    keys {} \
    flags {}

  if { $args == {} } {
    utl::error GPL 94 "placement_cluster requires a list of instances."
  }

  if { [llength $args] == 1 } {
    set args [lindex $args 0]
  }

  set insts []
  foreach inst_name $args {
    lappend insts {*}[gpl::parse_inst_names placement_cluster $inst_name]
  }
  utl::info GPL 96 "Created placement cluster of [llength $insts] instances."

  gpl::placement_cluster_cmd $insts
}

namespace eval gpl {
proc get_global_placement_uniform_density { args } {
  if { [ord::get_db_block] == "NULL" } {
    utl::error GPL 114 "No design block found."
  }

  sta::parse_key_args "get_global_placement_uniform_density" args \
    keys { -pad_left -pad_right } \
    flags {} ;# checker off

  # no need for init IP, TD and RD
  gpl::set_initial_place_max_iter_cmd 0
  gpl::set_routability_driven_mode 0
  gpl::set_timing_driven_mode 0


  # pad setting
  if { [info exists keys(-pad_left)] } {
    set pad_left $keys(-pad_left)
    sta::check_positive_integer "-pad_left" $pad_left
    gpl::set_pad_left_cmd $pad_left
  }
  if { [info exists keys(-pad_right)] } {
    set pad_right $keys(-pad_right)
    sta::check_positive_integer "-pad_right" $pad_right
    gpl::set_pad_right_cmd $pad_right
  }

  set uniform_density 0
  if { [ord::db_has_rows] } {
    sta::check_argc_eq0 "get_global_placement_uniform_density" $args

    set uniform_density [gpl::get_global_placement_uniform_density_cmd]
    gpl::replace_reset_cmd
  } else {
    utl::error GPL 131 "No rows defined in design. Use initialize_floorplan to add rows."
  }
  return $uniform_density
}

proc parse_inst_names { cmd names } {
  set dbBlock [ord::get_db_block]
  set inst_list {}
  foreach inst [get_cells $names] {
    lappend inst_list [sta::sta_to_db_inst $inst]
  }

  if { [llength $inst_list] == 0 } {
    utl::error GPL 95 "Instances {$names} for $cmd command were not found."
  }

  return $inst_list
}
}
