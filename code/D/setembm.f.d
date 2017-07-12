setembm.f
size.h
param.h
pconst.h
stdunits.h
calendar.h
solve.h
switch.h
coord.h
grdvar.h
cembm.h
atm.h
insolation.h
ice.h
evp.h
riv.h
tmngr.h
levind.h
csbc.h
scalar.h
veg.h
#if defined O_embm
# if defined O_ice
#  if defined O_ice_cpts
#  endif
#  if defined O_ice_evp
#  endif
# endif
# if defined O_embm_annual
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_embm_adiff
# endif
# if defined O_carbon_14_coupled
# endif
# if defined O_embm_explicit
# endif
# if defined O_time_averages
# endif
# if  defined O_co2emit_data || defined O_co2emit_data_transient
#  if defined O_co2emit_data_transient
#  endif
#  if defined O_carbon_co2_2d
#  endif
# endif
# if  defined O_co2ccn_data || defined O_co2ccn_data_transient
#  if defined O_co2ccn_data_transient
#  endif
# endif
# if !defined O_carbon_co2_2d
# endif
# if !defined O_orbit_user
#  if defined O_orbit_transient
#  endif
# endif
# if defined O_embm_explicit
# endif
# if defined O_embm_annual
# endif
# if !defined O_embm_explicit
#  if defined O_global_sums
#   if defined O_embm_slap
#   else
#   endif
#  else
#  endif
#  if defined O_embm_essl
#  endif
#  if defined O_embm_sparskit
#  endif
#  if defined O_embm_slap
#  endif
#  if defined O_embm_mgrid
#  endif
#  if defined O_embm_adi
#  endif
# else
# endif
# if defined O_ice
# else
# endif
# if defined O_embm_solve2x || defined O_embm_solve2y
#  if defined O_embm_solve2y
#  else
#  endif
#  if defined O_embm_solve2x
#  else
#  endif
#  if defined O_embm_solve2x
#  endif
#  if defined O_embm_solve2y
#  endif
#  if defined O_embm_solve2x && defined O_embm_solve2y
#  endif
#  if defined O_embm_solve2x
#  endif
#  if defined O_embm_solve2y
#  endif
#  if defined O_embm_solve2x && defined O_embm_solve2y
#  endif
# endif
# if defined O_embm_solve2y
# else
# endif
# if defined O_embm_solve2x
# else
# endif
# if defined O_ice_cpts && defined O_ice
#  if !defined O_ice_cpts5 && !defined O_ice_cpts10 && defined O_roth_press
#  endif
#  if defined O_ice_cpts3
#  elif defined O_ice_cpts5
#  elif defined O_ice_cpts10
#  endif
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_carbon && defined O_carbon_co2_2d
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
# endif
# if !defined O_mom
# endif
# if defined O_landice_data
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if !defined O_embm_explicit
#  if defined O_embm_slap
#  elif defined O_embm_essl
#  elif defined O_embm_sparskit
#  endif
# endif
# if defined O_landice_data
# else
# endif
# if defined O_sealev_data
# endif
# if defined O_restart_2
# endif
# if !defined O_mom
# endif
# if defined O_sealev
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_embm_awind || defined O_embm_adiff
#  if defined O_embm_awind || defined O_embm_adiff
#  else
#  endif
#  if defined O_embm_adiff
#  endif
# endif
#  if defined O_ice && !defined O_ice_cpts
#  endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_embm_solve2y
# endif
# if defined O_embm_solve2x
# endif
# if defined O_time_averages
# endif
# if defined O_time_step_monitor
# endif
#endif
