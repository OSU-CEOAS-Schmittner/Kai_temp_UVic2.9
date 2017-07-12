param.h
#if defined O_biharmonic
#endif
#if defined O_symmetry
#else
#endif
#if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
# if defined O_pressure_gradient_average
#  if defined O_biharmonic || defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker
#  else
#  endif
# else
# endif
#else
#endif
