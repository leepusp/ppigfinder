# ppigFinder v29_11_7 hardcoded inventory

Source: `legacy_sources/ppigfinder_v29_11_7.py`

Lines: `26508`

## Imports

- `PyQt5.QtCore`
- `PyQt5.QtGui`
- `PyQt5.QtPrintSupport`
- `PyQt5.QtSvg`
- `PyQt5.QtWidgets`
- `PyQt6.QtCore`
- `PyQt6.QtGui`
- `PyQt6.QtPrintSupport`
- `PyQt6.QtSvg`
- `PyQt6.QtWidgets`
- `base64`
- `collections`
- `copy`
- `csv`
- `datetime`
- `gc`
- `gzip`
- `hashlib`
- `ijson`
- `importlib`
- `importlib.metadata`
- `io`
- `json`
- `math`
- `matplotlib`
- `matplotlib.backends.backend_agg`
- `matplotlib.backends.backend_pdf`
- `matplotlib.backends.backend_qt5agg`
- `matplotlib.backends.backend_qtagg`
- `matplotlib.colors`
- `matplotlib.figure`
- `matplotlib.patches`
- `matplotlib.pyplot`
- `numpy`
- `os`
- `paramiko`
- `pathlib`
- `platform`
- `pyrodigal`
- `random`
- `re`
- `scipy.ndimage`
- `shlex`
- `shutil`
- `stat`
- `struct`
- `subprocess`
- `sys`
- `tempfile`
- `threading`
- `time`
- `traceback`
- `xml.etree`
- `xml.etree.ElementTree`

## Classes

### `AdvancedORFAnalyzer` line 1648

- `translate`
- `_translate_pyrodigal_gene`
- `reverse_complement`
- `gc_content`
- `find_orfs`
- `find_orfs_pyrodigal`
- `find_orfs_hybrid`
- `classify_domains`
- `analyze_neighborhood`
- `calc_evalue`
- `_blosum_score`
- `kmer_blast`
- `run_ncbi_blast`
- `_parse_blast_alignments`
- `sw_blast`
- `hmm_scan_orfs`
- `_resolve_pfam_executable`
- `_to_wsl_path`
- `pfam_prepare_library`
- `_tail`
- `_hmmscan_error_detail`
- `pfam_library_scan`
- `pfam_scan_with_alignments`
- `_hmmer3_search`
- `_parse_stockholm_aln`
- `_pssm_scan`
- `_parse_hmm_to_pssm`
- `generate_alphafold_input`
- `score_candidate`

### `AnalysisWorker` line 4054

- `__init__`
- `run`

### `GenomeMapWidget` line 4078

- `__init__`
- `set_data`
- `set_zoom`
- `_clamp_pan`
- `wheelEvent`
- `mousePressEvent`
- `mouseReleaseEvent`
- `mouseMoveEvent`
- `resizeEvent`
- `paintEvent`

### `_PpiAnalysisCanvas` line 4499

- `__init__`
- `_draw_axes`
- `export_svg`
- `export_png`
- `paintEvent`
- `_hit_test`
- `mousePressEvent`
- `mouseMoveEvent`

### `_PpiScatterWidget` line 4714

- `__init__`
- `set_data`
- `set_data_generic`
- `paint`

### `_PpiHistWidget` line 4858

- `__init__`
- `set_data`
- `paint`

### `_PpiRocWidget` line 4910

- `__init__`
- `set_data`
- `paint`

### `_Ppi2DScatterWidget` line 4961

- `__init__`
- `set_data`
- `set_data_generic`
- `paint`

### `_PpiNeighDensityWidget` line 5115

- `__init__`
- `set_data`
- `paint`

### `_PpiDomainHeatmapWidget` line 5150

- `__init__`
- `set_data`
- `paint`

### `_PpiArcMapWidget` line 5236

- `__init__`
- `set_data`
- `_genomic_to_x`
- `_x_to_genomic`
- `_arc_colour`
- `_orf_colour`
- `paintEvent`
- `_arc_hit`
- `_orf_hit`
- `mouseMoveEvent`
- `mousePressEvent`
- `mouseReleaseEvent`
- `wheelEvent`

### `_OrfNumericItem` line 5629

- `__init__`
- `_numeric_key`
- `__lt__`

### `WindowManager` line 5686

- `__init__`
- `_available_screen_rect`
- `apply_default_size`
- `apply_relative`
- `apply_fixed`
- `apply_maximized`
- `apply_fullscreen`
- `center_on_screen`
- `save_geometry`
- `restore_geometry`
- `register_subwindow`
- `open_custom_size_dialog`

### `CustomSizeDialog` line 5825

- `__init__`
- `_on_preset_changed`
- `_update_pixel_preview`
- `_apply_special`
- `get_result`

### `_DetachedTabWindow` line 5962

- `__init__`
- `closeEvent`

### `DetachableTabWidget` line 6060

- `__init__`
- `_on_tab_double_clicked`
- `mouseDoubleClickEvent`
- `addTab`
- `enable_body_double_click_detach`
- `eventFilter`
- `_on_tab_context_menu`
- `detach_tab`
- `_on_detached_closed`

### `AF3CommandBuilder` line 6315

- `detect_placeholders`
- `auto_context`
- `resolve`
- `_resolve_computed`
- `shell_quote_context`

### `AF3ProfileManager` line 6403

- `__init__`
- `load`
- `clear_user_profiles`
- `save`
- `get`
- `names`
- `add_or_update`
- `delete`
- `duplicate`

### `FlexibleAF3SubmitWidget` line 6616

- `__init__`
- `_build_ui`
- `_refresh_profile_combo`
- `_on_profile_changed`
- `_load_profile`
- `_build_param_widgets`
- `_on_template_changed`
- `_refresh_preview`
- `_on_new_profile`
- `_on_duplicate_profile`
- `_on_save_profile`
- `_on_delete_profile`
- `_on_import_profile`
- `_on_export_profile`
- `_on_reset_profiles`
- `_on_popout`
- `build_command`
- `current_profile_summary`
- `get_state`
- `set_state`

### `ppigFinderApp` line 7237

- `__init__`
- `_setup_ui`
- `_apply_accessible_ui_defaults`
- `_prepare_table`
- `_create_menus`
- `_build_window_menu`
- `_remove_panels_menu_action`
- `_build_panels_menu`
- `_core_right_panel_indices_for_top_buttons`
- `_compact_right_panel_label`
- `_build_top_right_panel_buttons`
- `_refresh_top_right_panel_button_state`
- `_panel_tooltip_for`
- `_apply_top_menu_tooltips`
- `_switch_right_view_by_index`
- `closeEvent`
- `_create_toolbar`
- `_create_statusbar`
- `_create_central`
- `_init_filter_vars`
- `_create_center_panel`
- `_create_right_navigation_bar`
- `_create_right_panel`
- `_right_view_label`
- `_sync_right_view_selector`
- `_on_right_view_combo_changed`
- `_on_right_tab_changed`
- `_switch_right_view_by_label`
- `_detach_current_right_view`
- `_redock_detached_right_views`
- `_create_genome_tab`
- `_create_blast_query_tab`
- `_create_blast_results_tab`
- `_create_dna_tab`
- `_create_protein_tab`
- `_create_domains_tab`
- `_create_neighborhood_tab`
- `_create_hmm_tab`
- `_update_hmm_profile_table`
- `_hmm_profile_anchor_name`
- `_on_hmm_profile_row_click`
- `_on_hmm_profile_right_click`
- `_af3_quick_homodimer_many`
- `_remove_hmm_profile`
- `_edit_hmm_profile`
- `_create_af3_tab`
- `_show_action_required`
- `_show_operation_error`
- `_run_worker`
- `analyze_orfs`
- `analyze_orfs_pyrodigal`
- `analyze_orfs_hybrid`
- `classify_all_domains`
- `annotate_pfam_library`
- `_apply_pfam_hits`
- `_pfam_search_orf_selection`
- `_pfam_search_selection`
- `_merge_pfam_sel_into_annotation`
- `_pfam_alignment_html`
- `_render_pfam_selection_results`
- `_archive_pfam_selection_to_tab`
- `_build_blast_cmd_string`
- `_refresh_blast_cmd_preview`
- `_run_custom_blast_cmd`
- `run_blast`
- `_show_blast_results`
- `_color_aln_html`
- `_color_mid_html`
- `_parse_fasta_query`
- `_validate_genome_sequence`
- `load_fasta`
- `load_multi_fasta`
- `load_snapgene`
- `load_genbank`
- `_on_sequence_loaded`
- `load_hmm`
- `_add_hmm_file`
- `save_fasta`
- `_orf_table_headers`
- `_orf_table_rows_tsv`
- `_orf_table_copy_selection`
- `_orf_table_copy_rows`
- `_export_orf_table`
- `_export_orf_fasta`
- `save_report_tsv`
- `_build_manifest`
- `save_project`
- `save_project_as`
- `load_project`
- `_script_run_directory`
- `export_map_pdf`
- `export_vector_figure`
- `_export_svg_figure`
- `_figure_saved_msg`
- `export_genbank`
- `export_snapgene`
- `_update_orfs_list`
- `_on_orf_note_changed`
- `_update_info`
- `_update_map`
- `_set_zoom`
- `_on_map_zoom_changed`
- `filter_orfs`
- `_on_text_orf_click`
- `_on_legend_click`
- `_update_hits_legend`
- `_on_orf_table_select`
- `_on_map_orf_click`
- `_parse_orf_idx_from_text`
- `_show_orf_details`
- `search_orf_in_map`
- `_select_and_center_orf`
- `_on_orf_right_click`
- `_switch_to_af3_tab`
- `_orf_table_resolve_sel_indices`
- `_orf_table_avg_residues`
- `_orf_table_af3_suggest_mode`
- `_orf_table_af3_add_multi`
- `_orf_table_af3_predict_allvsall`
- `_orf_table_af3_quick_neighbors`
- `_orf_table_af3_quick_homodimer`
- `_annotate_orf`
- `_color_orf`
- `_copy_to_clipboard`
- `_paste_clipboard`
- `_clear_query`
- `_validate_query`
- `load_blast_query_fasta`
- `_copy_blast_hit`
- `_copy_blast_all`
- `_save_blast_results`
- `_analyze_neighborhood`
- `_add_hmm_profile`
- `_add_hmm_multi`
- `_add_hmm_from_folder`
- `_hmm_search_all`
- `_show_hmm_results`
- `_af3_clean_seq`
- `_af3_validate_seq`
- `_af3_custom_dialog`
- `_af3_add_custom_sequence`
- `_af3_add_custom_row`
- `_af3_find_custom_entry`
- `_af3_edit_row_sequence`
- `_af3_resolve_entry_sequence`
- `_af3_seq_block`
- `_af3_add_selected`
- `_af3_add_orf_by_index`
- `_af3_add_hmm_hits`
- `_af3_add_all_orfs`
- `_af3_remove_orf`
- `_af3_clear_all`
- `_af3_show_job_size_config`
- `_af3_fragment_ranges`
- `_af3_name_with_fragment_ranges`
- `_af3_partition_pair_job`
- `_af3_build_fragment_jobs`
- `_af3_finalize_generated_jobs`
- `_af3_generate_jobs`
- `_af3_update_jobs_table`
- `_af3_export_json`
- `_af3_export_json_batch`
- `_af3_export_colabfold`
- `_af3_export_slurm_array`
- `_af3_show_ranking`
- `_af3_clear_jobs`
- `_af3_rebuild_custom_rows`
- `_af3_add_custom_job`
- `_af3_jobs_right_click`
- `_af3_sel_table_click`
- `_af3_jobs_table_click`
- `_af3_update_predict_pair_btn`
- `_af3_predict_selected_pair`
- `_af3_sel_table_right_click`
- `_af3_predict_single_vs_neighbors`
- `_af3_export_selected_jobs_json`
- `_af3_delete_selected_jobs`
- `keyPressEvent`
- `_show_genome_map_params`
- `_show_ppi_plot_params`
- `_load_af3_results_from_menu`
- `_export_ppi_arc_map_svg`
- `_export_ppi_arc_map_tsv`
- `_export_ppi_analysis_svg`
- `_export_ppi_analysis_png`
- `_export_af3a_tsv`
- `_export_af3a_motifs_tsv`
- `_export_af3a_anchor_windows`
- `_select_right_panel_by_label`
- `_show_orf_params`
- `_show_blast_params`
- `_show_hmm_params`
- `_show_manual`
- `_show_tutorial`
- `_show_install`
- `_show_system_check`
- `_show_help_dlg`
- `_show_help_interaction_results`
- `_show_help_ppi_map`
- `_show_help_ppi_analysis`
- `_show_help_references`
- `_show_about`
- `_af3a_user_role`
- `_create_af3_analysis_tab`
- `_af3a_load_folder`
- `_af3a_clear`
- `_af3a_classify_json`
- `_af3a_find_job_files`
- `_af3a_discover_job_dirs`
- `_af3a_scan_folder`
- `_af3a_seq_fingerprint`
- `_af3a_seq_seed_fingerprint`
- `_af3a_completeness_status`
- `_af3a_detect_truncation`
- `_af3a_assign_duplicate_groups`
- `_af3a_show_validation_report`
- `_af3a_parse_job`
- `_af3a_derive_token_plddt`
- `_af3a_chains_from_cif`
- `_af3a_label_cc`
- `_af3a_binary_closing`
- `_compute_pae_hotspot`
- `_af3a_visible_results`
- `_af3a_resolve_orf_ref`
- `_af3a_chain_domains`
- `_af3a_chain_offset_aa`
- `_af3a_pae_color`
- `_af3a_locate_anchor_window`
- `_af3a_anchor_windows_visible_results`
- `_af3a_build_anchor_window_table`
- `_af3a_export_anchor_windows`
- `_af3a_build_selected_interface_schematic`
- `_af3a_detect_motifs_one_pair`
- `_af3a_detect_motifs`
- `_af3a_contact_str`
- `_af3a_extract_orf_order_value`
- `_af3a_sort_by_orf_order`
- `_af3a_set_hotspot_radius`
- `_af3a_seq_map_for_result`
- `_af3a_anchor_sequence_texts`
- `_af3a_update_anchor_sequence_columns`
- `_af3a_populate_table`
- `_af3a_table_right_click`
- `_af3a_invert_selection`
- `_af3a_select_high_confidence`
- `_af3a_select_no_chain_info`
- `_af3a_select_duplicates`
- `_af3a_select_with_attribute`
- `_af3a_select_completeness`
- `_af3a_select_seq_status`
- `_af3a_remove_all_true_duplicates`
- `_af3a_res_for_view_row`
- `_af3a_copy_selected_job_names`
- `_af3a_copy_selected_rows_tsv`
- `_af3a_delete_selected_rows`
- `_af3a_unload_heavy_except`
- `_af3a_ensure_heavy_loaded`
- `_af3a_on_select`
- `_open_af3_job_folder`
- `_af3a_on_double_click`
- `_af3a_open_many`
- `_af3a_replot_selected`
- `_af3a_clear_plots`
- `_af3a_plot_job`
- `_af3a_canvas_context_menu`
- `_af3a_export_pdf`
- `_af3a_apply_filters`
- `_af3a_export_tsv`
- `_af3a_selected_result_index`
- `_af3a_refresh_hotspot_status`
- `_af3a_rerun_motifs`
- `_af3a_run_hotspots_selected`
- `_af3a_visible_ranked_candidates`
- `_af3a_run_hotspots_visible_top`
- `_af3a_run_hotspots_for_results`
- `_af3a_export_motifs_tsv`
- `_af3a_build_motif_table`
- `_af3a_highlight_motif`
- `_ppi_arc_update_legend`
- `_ppi_thresholds_changed`
- `_ppi_an_on_dot_click`
- `_ppi_predictions_for_orf`
- `_ppi_open_af3_prediction_for_orf`
- `_ppi_an_on_dot_context`
- `_ppi_arc_set_labels`
- `_ppi_arc_set_bidir`
- `_ppi_arc_set_topn`
- `_create_ppi_analysis_tab`
- `_ppi_analysis_current_widget`
- `_ppi_analysis_export_svg`
- `_ppi_analysis_export_png`
- `_ppi_detect_query`
- `_ppi_populate_posset`
- `_ppi_analysis_refresh`
- `_create_ppi_arc_map_tab`
- `_ppi_arc_map_refresh`
- `_ppi_arc_on_click`
- `_ppi_arc_on_orf_click`
- `_ppi_arc_on_hover`
- `_ppi_arc_map_export_svg`
- `_ppi_arc_map_export_tsv`
- `_create_hpc_server_tab`
- `_dv_build_connect_tab`
- `_dv_build_submit_tab`
- `_dv_build_monitor_tab`
- `_dv_cancel_selected_job`
- `_dv_export_job_history`
- `_dv_build_results_tab`
- `_dv_log`
- `_dv_set_connected`
- `_hpc_connect`
- `_dv_on_connected`
- `_dv_on_connect_error`
- `_dv_on_scheduler_changed`
- `_dv_on_env_method_changed`
- `_dv_build_activation_prefix`
- `_dv_scheduler_cmds`
- `_dv_load_module`
- `_dv_on_module_done`
- `_dv_on_module_error`
- `_hpc_disconnect`
- `_dv_ssh_exec`
- `_dv_remote_child_path_ok`
- `_dv_load_from_session`
- `_dv_show_batch_group_config`
- `_dv_load_from_session_batch`
- `_dv_load_from_files`
- `_dv_clear_submit_list`
- `_dv_download_staged_jsons`
- `_dv_af3_apply_preset`
- `_dv_refresh_cmd_preview`
- `_dv_upload_only`
- `_dv_upload_and_submit`
- `_dv_do_upload`
- `_dv_on_upload_done`
- `_dv_update_poll_timer`
- `_hpc_poll_queue`
- `_dv_on_poll_done`
- `_dv_refresh_monitor_table`
- `_dv_cancel_job_by_id`
- `_dv_autofill_results_path`
- `_dv_browse_remote_output`
- `_dv_choose_local_dest`
- `_dv_list_remote_output`
- `_dv_on_list_done`
- `_dv_download_results`
- `_dv_sftp_get_dir`
- `_dv_on_download_done`
- `_dv_import_ranking_from_path`

### `_NumericItem` line 20239

- `__init__`
- `__lt__`

### `_HE` line 17238

- `eventFilter`


## Top-level functions

- `detect_backends` line 992
- `_package_status` line 1083
- `build_system_check_report` line 1097
- `format_system_check_report` line 1163
- `_snapgene_read_packets` line 1221
- `parse_snapgene_dna` line 1240
- `_sg_packet` line 1329
- `write_snapgene_dna` line 1334
- `parse_genbank` line 1419
- `write_genbank` line 1559
- `t` line 3535
- `_setup_emoji_font` line 26339
- `_apply_text_fallback` line 26397
- `_check_dependencies_at_startup` line 26410
- `main` line 26475
- `_loc` line 1585
- `_wrap_qual` line 1592
- `translate` line 1719
- `_translate_pyrodigal_gene` line 1733
- `reverse_complement` line 1749
- `gc_content` line 1753
- `find_orfs` line 1757
- `find_orfs_pyrodigal` line 1795
- `find_orfs_hybrid` line 1900
- `classify_domains` line 1996
- `analyze_neighborhood` line 2006
- `calc_evalue` line 2054
- `_blosum_score` line 2063
- `kmer_blast` line 2070
- `run_ncbi_blast` line 2204
- `_parse_blast_alignments` line 2314
- `sw_blast` line 2439
- `hmm_scan_orfs` line 2528
- `_resolve_pfam_executable` line 2552
- `_to_wsl_path` line 2562
- `pfam_prepare_library` line 2568
- `_tail` line 2634
- `_hmmscan_error_detail` line 2638
- `pfam_library_scan` line 2671
- `pfam_scan_with_alignments` line 2811
- `_hmmer3_search` line 2979
- `_parse_stockholm_aln` line 3142
- `_pssm_scan` line 3210
- `_parse_hmm_to_pssm` line 3251
- `generate_alphafold_input` line 3293
- `score_candidate` line 3302
- `__init__` line 4060
- `run` line 4066
- `__init__` line 4083
- `set_data` line 4099
- `set_zoom` line 4113
- `_clamp_pan` line 4131
- `wheelEvent` line 4139
- `mousePressEvent` line 4155
- `mouseReleaseEvent` line 4168
- `mouseMoveEvent` line 4171
- `resizeEvent` line 4180
- `paintEvent` line 4186
- `__init__` line 4512
- `_draw_axes` line 4532
- `export_svg` line 4597
- `export_png` line 4614
- `paintEvent` line 4631
- `_hit_test` line 4638
- `mousePressEvent` line 4653
- `mouseMoveEvent` line 4683
- `__init__` line 4716
- `set_data` line 4726
- `set_data_generic` line 4741
- `paint` line 4752
- `__init__` line 4859
- `set_data` line 4868
- `paint` line 4877
- `__init__` line 4912
- `set_data` line 4921
- `paint` line 4936
- `__init__` line 4969
- `set_data` line 4978
- `set_data_generic` line 4987
- `paint` line 4997
- `__init__` line 5116
- `set_data` line 5124
- `paint` line 5129
- `__init__` line 5159
- `set_data` line 5166
- `paint` line 5176
- `__init__` line 5254
- `set_data` line 5287
- `_genomic_to_x` line 5297
- `_x_to_genomic` line 5308
- `_arc_colour` line 5320
- `_orf_colour` line 5332
- `paintEvent` line 5346
- `_arc_hit` line 5538
- `_orf_hit` line 5551
- `mouseMoveEvent` line 5559
- `mousePressEvent` line 5577
- `mouseReleaseEvent` line 5601
- `wheelEvent` line 5608
- `__init__` line 5638
- `_numeric_key` line 5643
- `__lt__` line 5661
- `__init__` line 5697
- `_available_screen_rect` line 5706
- `apply_default_size` line 5719
- `apply_relative` line 5732
- `apply_fixed` line 5747
- `apply_maximized` line 5755
- `apply_fullscreen` line 5758
- `center_on_screen` line 5761
- `save_geometry` line 5768
- `restore_geometry` line 5777
- `register_subwindow` line 5790
- `open_custom_size_dialog` line 5808
- `__init__` line 5828
- `_on_preset_changed` line 5912
- `_update_pixel_preview` line 5925
- `_apply_special` line 5931
- `get_result` line 5935
- `__init__` line 5970
- `closeEvent` line 6050
- `__init__` line 6067
- `_on_tab_double_clicked` line 6082
- `mouseDoubleClickEvent` line 6087
- `addTab` line 6102
- `enable_body_double_click_detach` line 6113
- `eventFilter` line 6127
- `_on_tab_context_menu` line 6145
- `detach_tab` line 6173
- `_on_detached_closed` line 6212
- `detect_placeholders` line 6325
- `auto_context` line 6336
- `resolve` line 6356
- `_resolve_computed` line 6372
- `shell_quote_context` line 6394
- `__init__` line 6406
- `load` line 6418
- `clear_user_profiles` line 6545
- `save` line 6558
- `get` line 6566
- `names` line 6572
- `add_or_update` line 6579
- `delete` line 6593
- `duplicate` line 6601
- `__init__` line 6625
- `_build_ui` line 6647
- `_refresh_profile_combo` line 6825
- `_on_profile_changed` line 6832
- `_load_profile` line 6839
- `_build_param_widgets` line 6870
- `_on_template_changed` line 6898
- `_refresh_preview` line 6919
- `_on_new_profile` line 6934
- `_on_duplicate_profile` line 6953
- `_on_save_profile` line 6973
- `_on_delete_profile` line 6990
- `_on_import_profile` line 7012
- `_on_export_profile` line 7034
- `_on_reset_profiles` line 7052
- `_on_popout` line 7083
- `build_command` line 7176
- `current_profile_summary` line 7187
- `get_state` line 7194
- `set_state` line 7203
- `__init__` line 7239
- `_setup_ui` line 7436
- `_apply_accessible_ui_defaults` line 7460
- `_prepare_table` line 7550
- `_create_menus` line 7565
- `_build_window_menu` line 7682
- `_remove_panels_menu_action` line 7759
- `_build_panels_menu` line 7805
- `_core_right_panel_indices_for_top_buttons` line 7809
- `_compact_right_panel_label` line 7841
- `_build_top_right_panel_buttons` line 7858
- `_refresh_top_right_panel_button_state` line 7926
- `_panel_tooltip_for` line 7938
- `_apply_top_menu_tooltips` line 7971
- `_switch_right_view_by_index` line 8051
- `closeEvent` line 8057
- `_create_toolbar` line 8071
- `_create_statusbar` line 8164
- `_create_central` line 8170
- `_init_filter_vars` line 8194
- `_create_center_panel` line 8213
- `_create_right_navigation_bar` line 8405
- `_create_right_panel` line 8451
- `_right_view_label` line 8520
- `_sync_right_view_selector` line 8536
- `_on_right_view_combo_changed` line 8556
- `_on_right_tab_changed` line 8563
- `_switch_right_view_by_label` line 8576
- `_detach_current_right_view` line 8585
- `_redock_detached_right_views` line 8590
- `_create_genome_tab` line 8600
- `_create_blast_query_tab` line 8634
- `_create_blast_results_tab` line 8711
- `_create_dna_tab` line 8740
- `_create_protein_tab` line 8754
- `_create_domains_tab` line 8768
- `_create_neighborhood_tab` line 8779
- `_create_hmm_tab` line 8801
- `_update_hmm_profile_table` line 8899
- `_hmm_profile_anchor_name` line 8922
- `_on_hmm_profile_row_click` line 8932
- `_on_hmm_profile_right_click` line 8978
- `_af3_quick_homodimer_many` line 9085
- `_remove_hmm_profile` line 9110
- `_edit_hmm_profile` line 9118
- `_create_af3_tab` line 9231
- `_show_action_required` line 9680
- `_show_operation_error` line 9689
- `_run_worker` line 9696
- `analyze_orfs` line 9710
- `analyze_orfs_pyrodigal` line 9740
- `analyze_orfs_hybrid` line 9805
- `classify_all_domains` line 9870
- `annotate_pfam_library` line 9915
- `_apply_pfam_hits` line 10066
- `_pfam_search_orf_selection` line 10124
- `_pfam_search_selection` line 10139
- `_merge_pfam_sel_into_annotation` line 10200
- `_pfam_alignment_html` line 10232
- `_render_pfam_selection_results` line 10302
- `_archive_pfam_selection_to_tab` line 10306
- `_build_blast_cmd_string` line 10343
- `_refresh_blast_cmd_preview` line 10362
- `_run_custom_blast_cmd` line 10366
- `run_blast` line 10369
- `_show_blast_results` line 10445
- `_color_aln_html` line 10550
- `_color_mid_html` line 10565
- `_parse_fasta_query` line 10577
- `_validate_genome_sequence` line 10596
- `load_fasta` line 10616
- `load_multi_fasta` line 10656
- `load_snapgene` line 10659
- `load_genbank` line 10670
- `_on_sequence_loaded` line 10680
- `load_hmm` line 10687
- `_add_hmm_file` line 10692
- `save_fasta` line 10700
- `_orf_table_headers` line 10726
- `_orf_table_rows_tsv` line 10732
- `_orf_table_copy_selection` line 10745
- `_orf_table_copy_rows` line 10757
- `_export_orf_table` line 10767
- `_export_orf_fasta` line 10873
- `save_report_tsv` line 10928
- `_build_manifest` line 10963
- `save_project` line 11155
- `save_project_as` line 11230
- `load_project` line 11434
- `_script_run_directory` line 11859
- `export_map_pdf` line 11877
- `export_vector_figure` line 11995
- `_export_svg_figure` line 12143
- `_figure_saved_msg` line 12168
- `export_genbank` line 12181
- `export_snapgene` line 12189
- `_update_orfs_list` line 12206
- `_on_orf_note_changed` line 12344
- `_update_info` line 12360
- `_update_map` line 12372
- `_set_zoom` line 12381
- `_on_map_zoom_changed` line 12386
- `filter_orfs` line 12390
- `_on_text_orf_click` line 12412
- `_on_legend_click` line 12422
- `_update_hits_legend` line 12437
- `_on_orf_table_select` line 12481
- `_on_map_orf_click` line 12496
- `_parse_orf_idx_from_text` line 12501
- `_show_orf_details` line 12518
- `search_orf_in_map` line 12533
- `_select_and_center_orf` line 12549
- `_on_orf_right_click` line 12609
- `_switch_to_af3_tab` line 12811
- `_orf_table_resolve_sel_indices` line 12816
- `_orf_table_avg_residues` line 12829
- `_orf_table_af3_suggest_mode` line 12837
- `_orf_table_af3_add_multi` line 12859
- `_orf_table_af3_predict_allvsall` line 12895
- `_orf_table_af3_quick_neighbors` line 12942
- `_orf_table_af3_quick_homodimer` line 12994
- `_annotate_orf` line 13025
- `_color_orf` line 13069
- `_copy_to_clipboard` line 13115
- `_paste_clipboard` line 13119
- `_clear_query` line 13123
- `_validate_query` line 13126
- `load_blast_query_fasta` line 13136
- `_copy_blast_hit` line 13142
- `_copy_blast_all` line 13146
- `_save_blast_results` line 13149
- `_analyze_neighborhood` line 13154
- `_add_hmm_profile` line 13236
- `_add_hmm_multi` line 13240
- `_add_hmm_from_folder` line 13248
- `_hmm_search_all` line 13264
- `_show_hmm_results` line 13299
- `_af3_clean_seq` line 13431
- `_af3_validate_seq` line 13444
- `_af3_custom_dialog` line 13471
- `_af3_add_custom_sequence` line 13551
- `_af3_add_custom_row` line 13566
- `_af3_find_custom_entry` line 13592
- `_af3_edit_row_sequence` line 13598
- `_af3_resolve_entry_sequence` line 13649
- `_af3_seq_block` line 13665
- `_af3_add_selected` line 13673
- `_af3_add_orf_by_index` line 13678
- `_af3_add_hmm_hits` line 13692
- `_af3_add_all_orfs` line 13714
- `_af3_remove_orf` line 13870
- `_af3_clear_all` line 13875
- `_af3_show_job_size_config` line 13881
- `_af3_fragment_ranges` line 13980
- `_af3_name_with_fragment_ranges` line 14002
- `_af3_partition_pair_job` line 14020
- `_af3_build_fragment_jobs` line 14097
- `_af3_finalize_generated_jobs` line 14133
- `_af3_generate_jobs` line 14174
- `_af3_update_jobs_table` line 14526
- `_af3_export_json` line 14572
- `_af3_export_json_batch` line 14593
- `_af3_export_colabfold` line 14629
- `_af3_export_slurm_array` line 14674
- `_af3_show_ranking` line 14920
- `_af3_clear_jobs` line 14958
- `_af3_rebuild_custom_rows` line 14967
- `_af3_add_custom_job` line 15024
- `_af3_jobs_right_click` line 15114
- `_af3_sel_table_click` line 15162
- `_af3_jobs_table_click` line 15172
- `_af3_update_predict_pair_btn` line 15182
- `_af3_predict_selected_pair` line 15194
- `_af3_sel_table_right_click` line 15269
- `_af3_predict_single_vs_neighbors` line 15338
- `_af3_export_selected_jobs_json` line 15388
- `_af3_delete_selected_jobs` line 15414
- `keyPressEvent` line 15424
- `_show_genome_map_params` line 15439
- `_show_ppi_plot_params` line 15516
- `_load_af3_results_from_menu` line 15702
- `_export_ppi_arc_map_svg` line 15718
- `_export_ppi_arc_map_tsv` line 15734
- `_export_ppi_analysis_svg` line 15749
- `_export_ppi_analysis_png` line 15763
- `_export_af3a_tsv` line 15777
- `_export_af3a_motifs_tsv` line 15791
- `_export_af3a_anchor_windows` line 15805
- `_select_right_panel_by_label` line 15820
- `_show_orf_params` line 15846
- `_show_blast_params` line 16049
- `_show_hmm_params` line 16135
- `_show_manual` line 16177
- `_show_tutorial` line 16180
- `_show_install` line 16183
- `_show_system_check` line 16186
- `_show_help_dlg` line 16237
- `_show_help_interaction_results` line 16256
- `_show_help_ppi_map` line 16447
- `_show_help_ppi_analysis` line 16557
- `_show_help_references` line 16674
- `_show_about` line 16930
- `_af3a_user_role` line 17003
- `_create_af3_analysis_tab` line 17010
- `_af3a_load_folder` line 17468
- `_af3a_clear` line 17476
- `_af3a_classify_json` line 17500
- `_af3a_find_job_files` line 17601
- `_af3a_discover_job_dirs` line 17743
- `_af3a_scan_folder` line 17796
- `_af3a_seq_fingerprint` line 17935
- `_af3a_seq_seed_fingerprint` line 17958
- `_af3a_completeness_status` line 17972
- `_af3a_detect_truncation` line 18000
- `_af3a_assign_duplicate_groups` line 18031
- `_af3a_show_validation_report` line 18098
- `_af3a_parse_job` line 18272
- `_af3a_derive_token_plddt` line 18986
- `_af3a_chains_from_cif` line 19077
- `_af3a_label_cc` line 19200
- `_af3a_binary_closing` line 19237
- `_compute_pae_hotspot` line 19265
- `_af3a_visible_results` line 19326
- `_af3a_resolve_orf_ref` line 19352
- `_af3a_chain_domains` line 19377
- `_af3a_chain_offset_aa` line 19399
- `_af3a_pae_color` line 19430
- `_af3a_locate_anchor_window` line 19444
- `_af3a_anchor_windows_visible_results` line 19533
- `_af3a_build_anchor_window_table` line 19548
- `_af3a_export_anchor_windows` line 19621
- `_af3a_build_selected_interface_schematic` line 19656
- `_af3a_detect_motifs_one_pair` line 19729
- `_af3a_detect_motifs` line 19927
- `_af3a_contact_str` line 20014
- `_af3a_extract_orf_order_value` line 20052
- `_af3a_sort_by_orf_order` line 20106
- `_af3a_set_hotspot_radius` line 20122
- `_af3a_seq_map_for_result` line 20139
- `_af3a_anchor_sequence_texts` line 20161
- `_af3a_update_anchor_sequence_columns` line 20199
- `_af3a_populate_table` line 20221
- `_af3a_table_right_click` line 20619
- `_af3a_invert_selection` line 20740
- `_af3a_select_high_confidence` line 20751
- `_af3a_select_no_chain_info` line 20769
- `_af3a_select_duplicates` line 20782
- `_af3a_select_with_attribute` line 20807
- `_af3a_select_completeness` line 20817
- `_af3a_select_seq_status` line 20827
- `_af3a_remove_all_true_duplicates` line 20837
- `_af3a_res_for_view_row` line 20891
- `_af3a_copy_selected_job_names` line 20911
- `_af3a_copy_selected_rows_tsv` line 20920
- `_af3a_delete_selected_rows` line 20937
- `_af3a_unload_heavy_except` line 21006
- `_af3a_ensure_heavy_loaded` line 21031
- `_af3a_on_select` line 21060
- `_open_af3_job_folder` line 21108
- `_af3a_on_double_click` line 21138
- `_af3a_open_many` line 21153
- `_af3a_replot_selected` line 21194
- `_af3a_clear_plots` line 21257
- `_af3a_plot_job` line 21275
- `_af3a_canvas_context_menu` line 21786
- `_af3a_export_pdf` line 21915
- `_af3a_apply_filters` line 21943
- `_af3a_export_tsv` line 21987
- `_af3a_selected_result_index` line 22137
- `_af3a_refresh_hotspot_status` line 22157
- `_af3a_rerun_motifs` line 22166
- `_af3a_run_hotspots_selected` line 22176
- `_af3a_visible_ranked_candidates` line 22192
- `_af3a_run_hotspots_visible_top` line 22220
- `_af3a_run_hotspots_for_results` line 22259
- `_af3a_export_motifs_tsv` line 22350
- `_af3a_build_motif_table` line 22427
- `_af3a_highlight_motif` line 22532
- `_ppi_arc_update_legend` line 22577
- `_ppi_thresholds_changed` line 22617
- `_ppi_an_on_dot_click` line 22663
- `_ppi_predictions_for_orf` line 22671
- `_ppi_open_af3_prediction_for_orf` line 22751
- `_ppi_an_on_dot_context` line 22827
- `_ppi_arc_set_labels` line 22894
- `_ppi_arc_set_bidir` line 22901
- `_ppi_arc_set_topn` line 22908
- `_create_ppi_analysis_tab` line 22916
- `_ppi_analysis_current_widget` line 23097
- `_ppi_analysis_export_svg` line 23107
- `_ppi_analysis_export_png` line 23122
- `_ppi_detect_query` line 23137
- `_ppi_populate_posset` line 23169
- `_ppi_analysis_refresh` line 23187
- `_create_ppi_arc_map_tab` line 23496
- `_ppi_arc_map_refresh` line 23716
- `_ppi_arc_on_click` line 23817
- `_ppi_arc_on_orf_click` line 23840
- `_ppi_arc_on_hover` line 23844
- `_ppi_arc_map_export_svg` line 23860
- `_ppi_arc_map_export_tsv` line 23883
- `_create_hpc_server_tab` line 23916
- `_dv_build_connect_tab` line 23970
- `_dv_build_submit_tab` line 24129
- `_dv_build_monitor_tab` line 24338
- `_dv_cancel_selected_job` line 24448
- `_dv_export_job_history` line 24489
- `_dv_build_results_tab` line 24506
- `_dv_log` line 24577
- `_dv_set_connected` line 24586
- `_hpc_connect` line 24612
- `_dv_on_connected` line 24644
- `_dv_on_connect_error` line 24662
- `_dv_on_scheduler_changed` line 24671
- `_dv_on_env_method_changed` line 24689
- `_dv_build_activation_prefix` line 24702
- `_dv_scheduler_cmds` line 24724
- `_dv_load_module` line 24757
- `_dv_on_module_done` line 24795
- `_dv_on_module_error` line 24859
- `_hpc_disconnect` line 24868
- `_dv_ssh_exec` line 24884
- `_dv_remote_child_path_ok` line 24893
- `_dv_load_from_session` line 24922
- `_dv_show_batch_group_config` line 24951
- `_dv_load_from_session_batch` line 25032
- `_dv_load_from_files` line 25138
- `_dv_clear_submit_list` line 25183
- `_dv_download_staged_jsons` line 25191
- `_dv_af3_apply_preset` line 25244
- `_dv_refresh_cmd_preview` line 25254
- `_dv_upload_only` line 25366
- `_dv_upload_and_submit` line 25370
- `_dv_do_upload` line 25374
- `_dv_on_upload_done` line 25654
- `_dv_update_poll_timer` line 25882
- `_hpc_poll_queue` line 25891
- `_dv_on_poll_done` line 25912
- `_dv_refresh_monitor_table` line 25970
- `_dv_cancel_job_by_id` line 26054
- `_dv_autofill_results_path` line 26078
- `_dv_browse_remote_output` line 26094
- `_dv_choose_local_dest` line 26108
- `_dv_list_remote_output` line 26114
- `_dv_on_list_done` line 26170
- `_dv_download_results` line 26183
- `_dv_sftp_get_dir` line 26258
- `_dv_on_download_done` line 26269
- `_dv_import_ranking_from_path` line 26290
- `_fmt_bp` line 4263
- `_draw_shape` line 4803
- `_draw_shape` line 5046
- `_sub` line 6365
- `_is_legacy` line 6519
- `_on_dlg_close` line 7154
- `_toggle_fs` line 7725
- `_update_preview` line 9161
- `_open_picker` line 9177
- `_update_mode_desc` line 9460
- `_update_job_preview` line 9474
- `work` line 9726
- `done` line 9729
- `work` line 9766
- `done` line 9789
- `work` line 9841
- `done` line 9856
- `_tick` line 10021
- `work` line 10034
- `cleanup` line 10040
- `done` line 10048
- `on_error` line 10055
- `work` line 10170
- `done` line 10175
- `on_error` line 10192
- `work` line 10407
- `done` line 10433
- `_bg_write` line 11180
- `_step` line 11281
- `_do_write` line 11389
- `work` line 13272
- `done` line 13290
- `_refresh` line 13513
- `_update_count` line 13771
- `_upd` line 14757
- `_parse_orf_idx` line 15027
- `_pick_col` line 15615
- `_restore` line 15651
- `_toggle_filter` line 15976
- `_upd_info` line 15998
- `_export` line 16217
- `_make_hover` line 17236
- `_idx` line 17627
- `_load_small_json` line 17640
- `_summary_score` line 17650
- `_prefer_by_index` line 17707
- `_dir_has_summary` line 17755
- `_draw_track` line 19690
- `_ranges` line 20016
- `_c` line 20257
- `_fmt` line 22026
- `_fmt` line 22374
- `_plddt_color` line 22472
- `_on_motif_select` line 22522
- `_name_to_idx` line 22686
- `_key` line 22738
- `_other_label` line 22783
- `_leg_line` line 23671
- `_leg_rect` line 23679
- `_do_cancel` line 24473
- `_do_connect` line 24629
- `_do_module` line 24777
- `norm_remote` line 24902
- `_ctx_for` line 25299
- `_build_batch_list_for_job` line 25429
- `_upload_one_partition` line 25467
- `_submit_one` line 25496
- `_do_upload_submit` line 25561
- `_do_poll` line 25901
- `_do` line 26064
- `_do_list` line 26126
- `_do_download` line 26212
- `_on_cancel` line 10006
- `_tok_info` line 14380
- `_find_orfs_by_seq` line 18579
- `_xpos` line 19696
- `_resnum` line 19888
- `__init__` line 20240
- `__lt__` line 20245
- `_fmt_chain` line 20272
- `_on_hover` line 21623
- `_pae_right_click` line 21700
- `_plddt_right_click` line 21765
- `_arc_partner_indices` line 23274
- `_apply_transform` line 23304
- `_base_value` line 23326
- `_x_label_format` line 23427
- `_find_orf` line 23748
- `to_wsl_path` line 3007
- `_open_multi` line 12710
- `eventFilter` line 17239
- `_do_delete_resubmit` line 25796
- `_on_resubmit` line 25828
- `_chain_of` line 21635
- `_chain_letter` line 18490
- `_get` line 19123

## Keyword hits

### AlphaFold

16, 44, 137, 140, 167, 181, 187, 198, 297, 487, 502, 509, 568, 623, 640, 976, 3293, 3336, 3392, 3560, 3569, 3607, 3633, 3691, 3742, 3844, 3846, 3869, 6274, 6275, 6278, 6284, 6286, 6287, 6293, 6294, 6300, 6307, 6308, 6309, 6631, 7299, 7954, 7955, 7956, 7957, 8090, 8435, 8436, 8492, 8983, 9029, 9042, 9392, 11046, 11136, 11141, 11324, 12674, 12680, 12812, 12940, 12992, 13023, 13483, 13718, 14580, 14594, 14613, 14818, 15401, 16272, 16433, 16733, 16744, 16749, 16754, 16757, 16773, 16776

### AF3

16, 44, 53, 230, 232, 251, 298, 302, 321, 337, 340, 341, 342, 343, 344, 346, 347, 348, 352, 360, 363, 369, 370, 371, 372, 423, 476, 492, 501, 507, 510, 518, 519, 568, 573, 589, 610, 640, 649, 650, 664, 695, 706, 719, 721, 754, 755, 773, 849, 1074, 1078, 1079, 3336, 3392, 3425, 3426, 3427, 3429, 3430, 3431, 3432, 3433, 3434, 3435, 3440, 3441, 3442, 3443, 3444, 3445, 3446, 3447, 3448, 3449, 3450, 3451, 3452, 3453, 3454, 3455

### BLAST

15, 24, 119, 122, 157, 567, 995, 1015, 1032, 1111, 1137, 1139, 1651, 2070, 2071, 2202, 2204, 2205, 2206, 2217, 2236, 2239, 2254, 2287, 2295, 2296, 2298, 2303, 2306, 2309, 2314, 2315, 2319, 2335, 2361, 2378, 2435, 2439, 3343, 3379, 3383, 3389, 3390, 3414, 3415, 3416, 3417, 3418, 3419, 3420, 3421, 3465, 3466, 3467, 3468, 3469, 3470, 3471, 3472, 3473, 3474, 3475, 3476, 3493, 3494, 3495, 3496, 3497, 3498, 3499, 3500, 3501, 3502, 3503, 3504, 3505, 3522, 3523, 3524, 3525

### HMM

15, 26, 27, 28, 29, 30, 33, 39, 78, 126, 127, 129, 161, 293, 294, 307, 308, 311, 387, 397, 404, 416, 417, 448, 530, 536, 537, 540, 541, 546, 547, 548, 549, 552, 556, 557, 558, 559, 562, 563, 564, 568, 696, 753, 767, 995, 1015, 1112, 1140, 1142, 1144, 1624, 1652, 1665, 1666, 2526, 2528, 2529, 2531, 2532, 2533, 2537, 2540, 2541, 2544, 2547, 2548, 2549, 2553, 2555, 2556, 2558, 2569, 2571, 2572, 2573, 2577, 2579, 2582, 2583

### HMMER

26, 126, 127, 552, 995, 1015, 1112, 1140, 1142, 1144, 2529, 2531, 2532, 2553, 2555, 2556, 2558, 2577, 2609, 2610, 2629, 2666, 2689, 2740, 2741, 2838, 2893, 2979, 2980, 2985, 3134, 3211, 3252, 3256, 3263, 3283, 3365, 3513, 3515, 3658, 3672, 3792, 3793, 3814, 3901, 3911, 3920, 3929, 3935, 3936, 3937, 3943, 3944, 3957, 3959, 3971, 3976, 3977, 3983, 3996, 4006, 4011, 4029, 6297, 8011, 8019, 8033, 8156, 8158, 8615, 8624, 9927, 9929, 9931, 9932, 9935, 9936, 10144, 10146, 10147

### Pyrodigal

22, 108, 111, 145, 243, 737, 741, 754, 933, 934, 936, 1063, 1064, 1065, 1066, 1067, 1076, 1113, 1146, 1148, 1733, 1734, 1795, 1798, 1800, 1816, 1818, 1820, 1822, 1830, 1835, 1846, 1863, 1867, 1891, 1892, 1903, 1908, 1913, 1915, 1917, 1927, 1933, 1934, 1943, 1953, 3354, 3357, 3358, 3517, 3518, 3519, 3520, 3521, 3565, 3584, 3643, 3644, 3645, 3794, 3814, 3820, 3900, 3910, 3925, 3962, 3965, 3971, 3984, 3994, 4007, 4013, 4024, 4025, 4374, 4412, 7307, 8015, 8033, 8099

### Prodigal

112, 114, 146, 151, 1798, 3358, 3517, 15900, 16690, 16694, 16698

### SnapGene

13, 90, 250, 251, 710, 711, 1210, 1218, 1221, 1222, 1227, 1240, 1242, 1248, 1261, 1326, 1330, 1334, 1338, 1371, 1423, 1432, 1565, 1605, 1610, 3329, 3338, 3401, 7327, 7328, 7329, 10659, 10660, 10661, 10663, 10666, 10667, 10677, 11148, 11149, 11657, 11658, 11659, 11660, 12186, 12189, 12191, 12198, 12200

### GenBank

13, 90, 250, 251, 711, 1210, 1416, 1419, 1421, 1556, 1559, 1563, 1586, 3330, 3339, 3401, 10670, 10671, 10672, 10674, 12181, 12183, 12185, 12187

### QThread

872, 906, 4054

### QTimer

872, 906, 7336, 7392, 7394, 10019

### icon

5970, 5985, 5987, 5993, 6181, 6201, 6228, 7449, 8523, 8524, 20285, 20288, 20305, 20334, 23936, 23937, 23938, 24078, 24589, 24599, 26461, 26462

### toolbar

42, 64, 334, 336, 337, 339, 344, 354, 434, 536, 561, 623, 733, 867, 901, 3351, 3581, 3661, 3796, 3811, 3819, 3829, 4004, 7425, 7455, 7885, 7958, 8069, 8071, 8072, 8073, 8075, 8356, 9338, 15822, 16283, 16326, 16366, 16395, 16505, 16517, 16598, 16612, 16625, 16636, 16650, 16663, 16886, 16906, 16951, 17026, 17120, 17121, 18869, 18883, 19187, 22923, 23003, 23503

### backend

948, 951, 960, 989, 992, 999, 1035, 1038, 1049, 1057, 1059, 1061, 1064, 1073, 1111, 1112, 1113, 1114, 1137, 1140, 1146, 2206, 2531, 2558, 2577, 2689, 2838, 2985, 3511, 3788, 3795, 3796, 3813, 4009, 8154, 8155, 8156, 8157, 8204, 8613, 8614, 8615, 8616, 8622, 8623, 8624, 8625, 8626, 9929, 10144, 10235, 10358, 10411, 13301, 13302, 16061, 16931, 16932, 16933, 16968, 21283, 21932

### export

41, 42, 43, 44, 62, 83, 90, 346, 349, 350, 354, 355, 388, 409, 501, 502, 512, 524, 625, 678, 687, 711, 730, 733, 734, 735, 736, 762, 769, 1074, 1576, 3338, 3339, 3340, 3366, 3404, 3432, 3451, 3452, 3453, 3454, 3455, 3590, 3591, 3596, 3597, 3709, 3739, 3740, 3743, 3745, 3763, 3859, 3877, 4034, 4036, 4597, 4614, 6686, 6706, 6707, 6708, 7034, 7041, 7048, 7050, 7320, 7597, 7598, 7614, 7615, 7618, 7622, 7624, 7627, 7629, 7635, 7637, 7641, 7643

### HTML

311, 8735, 8890, 8923, 8943, 8944, 10232, 10233, 10234, 10236, 10237, 10238, 10240, 10242, 10244, 10245, 10254, 10262, 10267, 10272, 10282, 10289, 10290, 10293, 10296, 10298, 10299, 10300, 10304, 10329, 10446, 10449, 10452, 10453, 10454, 10455, 10456, 10458, 10461, 10462, 10463, 10465, 10473, 10479, 10480, 10481, 10482, 10496, 10497, 10499, 10503, 10530, 10531, 10532, 10535, 10536, 10539, 10540, 10541, 10546, 10547, 10550, 10565, 11023, 11026, 11146, 11171, 11172, 11174, 11382, 11384, 11666, 11667, 13300, 13305, 13306, 13307, 13308, 13309, 13311

### PDF

41, 42, 87, 90, 251, 350, 353, 354, 356, 388, 403, 409, 411, 464, 472, 3340, 3366, 3404, 3590, 3591, 3615, 3618, 3889, 7617, 7618, 7654, 7746, 7752, 7757, 8005, 8134, 8135, 11860, 11877, 11878, 11883, 11885, 11886, 11890, 11892, 11906, 11907, 11909, 11936, 11965, 12003, 12022, 12030, 12031, 12039, 12052, 12053, 12060, 12077, 12090, 12102, 16949, 16951, 16952, 16962, 17074, 17075, 17076, 17077, 17461, 21788, 21800, 21801, 21802, 21823, 21824, 21858, 21860, 21861, 21868, 21912, 21913, 21915, 21916, 21918

### SVG

83, 351, 4597, 4599, 4602, 4604, 7634, 7635, 7640, 7641, 11886, 11967, 11970, 11972, 11973, 11980, 12003, 12032, 12101, 12102, 12103, 12104, 12105, 12109, 12110, 12111, 12143, 12144, 12147, 12149, 12150, 12165, 12166, 15718, 15719, 15721, 15730, 15749, 15750, 15751, 15760, 16533, 16534, 16595, 16656, 21788, 21804, 21805, 21806, 21870, 21872, 21873, 21880, 22987, 22988, 22989, 22990, 23102, 23107, 23112, 23113, 23117, 23561, 23562, 23563, 23564, 23860, 23861, 23863, 23866, 23868, 23869, 23872, 23875

### SLURM

18, 51, 520, 728, 761, 762, 769, 772, 775, 3454, 3737, 3745, 6252, 6279, 6290, 7334, 7955, 9564, 9565, 13905, 14674, 14675, 14691, 14695, 14702, 14724, 14783, 14806, 14840, 14846, 14847, 14851, 14868, 14897, 14915, 14918, 24043, 24050, 24176, 24385, 24395, 24403, 24422, 24449, 24459, 24460, 24462, 24464, 24468, 24729, 24755, 24968, 25022, 25038, 25543, 25605, 25609, 25620, 25670, 25671, 25676, 25681, 25687, 25691, 25698, 25699, 25702, 25711, 25722, 25819, 25824, 25830, 25834, 25839, 25848, 25869, 25878, 25893, 25935, 25943

### sbatch

3745, 6252, 9566, 14677, 14687, 14703, 14829, 14830, 14831, 14832, 14833, 14834, 14835, 14836, 14837, 14842, 14880, 14898, 14905, 14913, 24030, 24031

### ssh

18, 51, 134, 772, 1077, 3999, 7332, 23927, 23976, 23989, 24023, 24114, 24450, 24474, 24575, 24613, 24630, 24646, 24654, 24760, 24786, 24860, 24869, 24874, 24875, 24878, 24884, 24886, 24888, 25387, 25407, 25481, 25541, 25565, 25571, 25576, 25798, 25815, 25888, 25894, 25902, 26056, 26065, 26116, 26128, 26132, 26185, 26214, 26218

### paramiko

133, 135, 243, 939, 940, 942, 1077, 3901, 3910, 3925, 3962, 3965, 3972, 3984, 3999, 23924, 23925, 23927, 23928, 24011, 24602, 24614, 24616, 24630, 24631, 24666, 26419

