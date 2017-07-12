fdifm.h
#if defined O_linearized_advection
#else
#endif
#if defined O_implicitvmix
#endif
#if defined O_consthmix && !defined O_biharmonic
# if defined O_neptune
# else
#  if defined O_anisotropic_viscosity
#  else
#  endif
# endif
#else
#endif
#if defined O_consthmix
# if defined O_biharmonic
# else
#  if defined O_neptune
#  else
#  endif
# endif
#else
# if defined O_smagnlmix
# endif
#endif
#if defined O_stream_function
# if defined O_damp_inertial_oscillation
# else
# endif
#else
#endif
