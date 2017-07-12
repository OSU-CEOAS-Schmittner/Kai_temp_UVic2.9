diag.h
#if defined O_mom
# if defined O_time_step_monitor
#  if defined O_tai_otsf
#  endif
#  if defined O_tai_slh
#  endif
# endif
# if defined O_energy_analysis
# endif
# if defined O_gyre_components
#  if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
#  endif
#  if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
#  else
#  endif
# endif
# if defined O_meridional_overturning
# endif
# if defined O_show_zonal_mean_of_sbc
# endif
# if defined O_tracer_yz
# endif
# if defined O_term_balances
# endif
#endif
