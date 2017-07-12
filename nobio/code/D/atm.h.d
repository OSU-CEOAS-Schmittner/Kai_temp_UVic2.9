atm.h
#if defined O_ice_evp || defined O_embm_awind
#else
#endif
#if defined O_embm_awind || defined O_embm_adiff
#endif
#if defined O_embm_awind
#endif
#if defined O_sulphate_data || defined O_sulphate_data_transient
#endif
#if defined O_carbon_co2_2d
# if defined O_co2emit_data || O_co2emit_data_transient
# endif
#endif
#if defined O_sealev || defined O_sealev_data
#endif
#if defined O_save_flxadj
#endif
#if defined O_co2emit_track_co2 || defined O_co2emit_track_co2_transient
#endif
#if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient || defined O_embm_vcs
#endif
#if defined O_time_averages
# if defined O_embm_awind
# endif
# if defined O_save_flxadj
# endif
# if defined O_save_embm_diff
# endif
# if defined O_carbon_co2_2d
# endif
# if defined O_carbon_co2_2d
#  if defined O_co2emit_data || defined O_co2emit_data_transient
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
#endif
