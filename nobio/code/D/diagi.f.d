diagi.f
size.h
param.h
pconst.h
stdunits.h
stab.h
ctmb.h
coord.h
ctavg.h
diag.h
iounit.h
switch.h
tmngr.h
#if defined O_mom
# if defined O_tracer_averages
# endif
# if defined O_stability_tests
# endif
# if defined O_energy_analysis
# endif
# if defined O_term_balances
# endif
# if defined O_gyre_components
# if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
# else
# endif
# endif
# if defined O_meridional_overturning
# endif
# if defined O_tracer_averages
# endif
# if defined O_time_step_monitor
# endif
# if defined O_show_zonal_mean_of_sbc
# endif
#endif
