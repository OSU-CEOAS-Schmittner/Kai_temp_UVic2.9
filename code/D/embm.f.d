embm.f
size.h
csbc.h
cembm.h
param.h
pconst.h
stdunits.h
atm.h
ice.h
#if defined O_embm
# if !defined O_mom
# endif
# if !defined O_mom
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_even_fluxes && defined O_mom
# endif
# if defined O_ice
# endif
# if defined O_time_averages
# endif
# if defined O_time_step_monitor
# endif
# if !defined O_mom
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_global_sums
# endif
# if defined O_global_sums
# endif
# if defined O_global_sums
# endif
# if defined O_mtlm
#  if defined O_landice_data
#  endif
#  if defined O_sealev || defined O_sealev_data
#  endif
# endif
# if defined O_ice_evp && defined O_ice
# elif defined O_embm_awind
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_global_sums
# endif
#endif
