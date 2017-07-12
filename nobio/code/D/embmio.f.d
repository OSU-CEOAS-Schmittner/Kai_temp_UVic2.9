embmio.f
size.h
param.h
pconst.h
stdunits.h
calendar.h
csbc.h
atm.h
solve.h
coord.h
grdvar.h
levind.h
ice.h
evp.h
mtlm.h
cembm.h
iounit.h
scalar.h
switch.h
tmngr.h
riv.h
cregin.h
#if defined O_embm
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_mtlm
# endif
# if defined O_orbit_transient
# endif
# if defined O_mom
# else
# endif
# if defined O_mtlm
# endif
# if defined O_embm_adiff
# endif
# if !defined O_mom
# endif
# if defined O_time_step_monitor
#  if defined O_mtlm
#  endif
#  if defined O_save_carbon_totals
#  else
#  endif
#  if !defined O_save_time_relyear0
#  endif
#  if defined O_save_time_endper
#  elif defined O_save_time_startper
#  else
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  else
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  endif
# endif
# if defined O_time_averages
#  if defined O_mtlm
#  endif
#  if defined O_landice_data
#  endif
#  if defined O_sealev || defined O_sealev_data
#  endif
#  if defined O_mtlm
#   if defined O_ice
#   endif
#  endif
#  if defined O_ice
#  endif
#  if !defined O_save_time_relyear0
#  endif
#  if defined O_save_time_endper
#  elif defined O_save_time_startper
#  else
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  else
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  endif
#  if defined O_save_embm_wind
#  endif
#  if defined O_embm_awind
#  endif
#  if defined O_mtlm
#  else
#  endif
#  if defined O_ice
#  endif
#  if defined O_ice_cpts && defined O_ice
#  endif
#  if defined O_ice_evp && defined O_ice
#  endif
#  if defined O_save_flxadj
#  endif
#  if defined O_save_embm_diff
#  endif
#  if defined O_landice_data
#  endif
#  if defined O_carbon_co2_2d
#   if defined O_co2emit_data || defined O_co2emit_data_transient
#   endif
#  endif
#  if defined O_sulphate_data || defined O_sulphate_data_transient
#  endif
# endif
# if defined O_orbit_transient
# endif
# if defined O_landice_data_transient
#  if defined O_ice && !defined O_ice_cpts
#  endif
#  if defined O_mtlm
#   if defined O_landice_data
#   else
#   endif
#  endif
#  if defined O_sealev
#  endif
# endif
# if defined O_embm_awind || defined O_embm_adiff
#  if defined O_embm_awind
#  endif
#  if defined O_embm_adiff
#  endif
# endif
#endif
#if defined O_embm && defined O_time_averages
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_save_embm_wind
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
#  if defined O_ice_cpts && defined O_ice
#  endif
#  if defined O_ice_evp && defined O_ice
#  endif
# endif
# if defined O_save_flxadj
# endif
# if defined O_save_embm_diff
# endif
# if defined O_landice_data
# endif
# if defined O_carbon_co2_2d
#  if defined O_co2emit_data || defined O_co2emit_data_transient
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_save_embm_wind
#  if defined O_carbon_co2_2d
#  endif
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
#  if defined O_ice_cpts
#  else
#  endif
#  if defined O_ice_evp && defined O_ice
#  endif
# endif
# if defined O_save_flxadj
# endif
# if defined O_save_embm_diff
# endif
# if defined O_landice_data
# endif
# if defined O_carbon_co2_2d
#  if defined O_co2emit_data || defined O_co2emit_data_transient
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_save_embm_wind
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
#  if defined O_ice_cpts
#  endif
#  if defined O_ice_evp && defined O_ice
#  endif
# endif
# if defined O_save_flxadj
# endif
# if defined O_save_embm_diff
# endif
# if defined O_landice_data
# endif
# if defined O_carbon_co2_2d
#  if defined O_co2emit_data || defined O_co2emit_data_transient
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
#endif
#if defined O_embm && defined O_time_step_monitor
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_tai_rad
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_embm_explicit
# else
# endif
# if defined O_landice_data
# endif
# if defined O_carbon_co2_2d
# else
# endif
# if defined O_carbon && defined O_carbon_13
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_tai_lo
# endif
# if defined O_ice
#  if defined O_ice_cpts
#  else
#  endif
# endif
# if defined O_tai_lo
# endif
# if defined O_ice
# endif
# if defined O_landice_data
# endif
# if defined O_tai_ns
#  if defined O_ice_cpts && defined O_ice
#  elif defined O_ice
#  endif
#  if defined O_ice
#  endif
#  if defined O_landice_data
#  endif
#  if defined O_ice_cpts && defined O_ice
#  elif defined O_ice
#  endif
#  if defined O_ice
#  endif
#  if defined O_landice_data
#  endif
# endif
# if defined O_tai_rad
#  if defined O_tai_lo
#  endif
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#  if defined O_carbon_13
#  endif
#  if defined O_npzd_iron
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
# endif
#endif
