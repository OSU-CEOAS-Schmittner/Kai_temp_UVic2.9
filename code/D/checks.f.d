checks.f
size.h
param.h
pconst.h
stdunits.h
accel.h
coord.h
csbc.h
grdvar.h
hmixc.h
iounit.h
levind.h
isopyc.h
mw.h
scalar.h
switch.h
vmixc.h
# if defined O_mom
# if defined O_isopycmix
# endif
# if defined O_xbts
# endif
# if defined O_pressure_gradient_average
# endif
# if defined O_linearized_advection
#  if defined O_fct
#  endif
#  if defined O_fourth_order_tracer_advection
#  endif
#  if defined O_quicker
#  endif
#  if !defined O_linearized_density
#  endif
# endif
# if defined O_consthmix
# endif
# if defined O_constvmix
# endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_bryan_lewis_horizontal
#  if !defined O_consthmix
#  endif
# endif
# if defined O_rigid_lid_surface_pressure
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  if defined O_hypergrid || defined O_oldrelax
#  endif
#  if defined O_sf_5_point
#  endif
# endif
# if !defined O_stream_function && !defined O_implicit_free_surface
#  if !defined O_rigid_lid_surface_pressure
#  endif
# endif
# if defined O_stream_function
#  if !defined O_sf_5_point && !defined O_sf_9_point
#  endif
#  if defined O_sf_5_point && defined O_sf_9_point
#  endif
# endif
# if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
# else
# endif
# if defined O_restorst
# else
# endif
# if defined O_fourth_order_tracer_advection
# endif
# if defined O_isopycmix
#  if defined O_consthmix && !defined O_biharmonic
#  endif
#  if defined O_biharmonic
#  endif
# else
#  if defined O_tidal_kv
#  endif
#  if defined O_gent_mcwilliams
#  endif
# endif
# if defined O_implicit_free_surface && defined O_stream_function
# endif
# if defined O_rigid_lid_surface_pressure && defined O_stream_function
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  if defined O_diagnostic_surf_height
#  endif
# endif
# if defined O_stream_function
#  if !defined O_conjugate_gradient && !defined O_oldrelax
#   if !defined O_hypergrid
#   endif
#  endif
#  if defined O_oldrelax && defined O_hypergrid
#  endif
#  if defined O_oldrelax && defined O_conjugate_gradient
#  endif
#  if defined O_oldrelax && defined O_conjugate_gradient
#  endif
#  if defined O_hypergrid && defined O_conjugate_gradient
#  endif
#  if defined O_hypergrid && defined O_conjugate_gradient
#  endif
#  if defined O_sf_9_point
#  endif
#  if defined O_sf_5_point
#  endif
#  if defined O_sf_9_point && defined O_oldrelax
#  endif
# endif
# if defined O_biharmonic && !defined O_consthmix
# endif
# if defined O_isopycmix
#  if defined O_consthmix
#  endif
# endif
# if defined O_meridional_tracer_budget
# else
# endif
# if defined O_xbts
# else
# endif
# if defined O_diagnostic_surf_height
# else
# endif
# if defined O_consthmix
#  if defined O_isopycmix
#  else
#  endif
#  if defined O_isopycmix
#  else
#  endif
# endif
# if defined O_shortwave
#  if !defined O_source_term
#  endif
#  if !defined O_embm
#  endif
# else
# endif
# if defined O_constvmix
# endif
# if defined O_implicitvmix
#  if defined O_fullconvect
#  endif
# else
#  if defined O_fullconvect
#  else
#  endif
# endif
# if !defined O_implicitvmix && !defined O_isopycmix
# endif
# if defined O_damp_inertial_oscillation
#  if defined O_fct
#  endif
# else
# endif
# if defined O_consthmix
# endif
# if defined O_cyclic && defined O_solid_walls
# endif
# if !defined O_cyclic && !defined O_solid_walls
# endif
# if defined O_solid_walls
# endif
# if !defined O_symmetry
# endif
# if !defined O_quicker && defined O_ncar_upwind3
# endif
# if defined O_fct
#  if defined O_fourth_order_tracer_advection
#  endif
#  if  defined O_quicker
#  endif
#  if defined O_fct_dlm1 && defined O_fct_dlm2
#  endif
#  if !defined O_fct_dlm1 && !defined O_fct_dlm2
#  endif
# else
#  if defined O_fct_dlm1 || defined O_fct_dlm2 || defined O_fct_3d
#  endif
# endif
#endif
