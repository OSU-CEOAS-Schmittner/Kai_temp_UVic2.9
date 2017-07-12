bardiv.f
#if defined O_mom
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  if defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  else
#  endif
#  if !defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  endif
#  if defined O_implicit_free_surface
#  else
#  endif
#  if defined O_remove_ps_checkerboard
#   if !defined O_implicit_free_surface
#   endif
#  endif
# endif
#endif
