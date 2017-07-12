clinic.f
size.h
param.h
pconst.h
stdunits.h
coord.h
csbc.h
grdvar.h
hmixc.h
emode.h
levind.h
mw.h
scalar.h
switch.h
vmixc.h
fdifm.h
diag.h
diaga.h
#if defined O_mom
# if defined O_neptune
# endif
# if defined O_biharmonic
# else
# endif
# if defined O_anisotropic_viscosity
# else
#  if defined O_consthmix
#  endif
# endif
# if defined O_pressure_gradient_average
# endif
# if defined O_pressure_gradient_average
# else
# endif
# if defined O_pressure_gradient_average
# else
# endif
# if !defined O_linearized_advection
#  if defined O_consthmix && !defined O_biharmonic
#   if defined O_anisotropic_viscosity
#   else
#   endif
#   if defined O_neptune
#   endif
#  endif
# endif
# if defined O_consthmix && defined O_biharmonic
# endif
# if defined O_smagnlmix && !defined O_consthmix
# endif
# if !defined O_linearized_advection
# endif
# if defined O_source_term || defined O_npzd || defined O_carbon_14
# endif
# if !defined O_linearized_advection
# endif
# if defined O_source_term || defined O_npzd || defined O_carbon_14
# endif
# if defined O_implicitvmix
# endif
# if defined O_symmetry
# endif
# if defined O_damp_inertial_oscillation
# else
# endif
# if defined O_fourfil || defined O_firfil
# endif
# if defined O_ice_evp
# endif
# if defined O_save_mixing_coeff
#  if !defined O_consthmix || defined O_biharmonic
#  else
#  endif
# endif
# if defined O_time_step_monitor
#  if defined O_symmetry
#  endif
# endif
# if defined O_energy_analysis
# endif
# if defined O_term_balances
# endif
# if defined O_xbts
# endif
# if defined O_energy_analysis
# endif
# if defined O_term_balances
# endif
# if defined O_xbts
# endif
#endif
#if defined O_mom && defined O_ice_evp
#endif
#if defined O_implicitvmix
# if defined O_xbts || defined O_energy_analysis || defined O_term_balances
# endif
# if defined O_xbts || defined O_energy_analysis || defined O_term_balances
# endif
#endif
