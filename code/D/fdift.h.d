fdift.h
#if defined O_linearized_advection
#else
# if defined O_fourth_order_tracer_advection || defined O_quicker
# else
#  if defined O_fct
#  else
#  endif
# endif
#endif
#if defined O_gent_mcwilliams && defined O_isopycmix
#endif
#if defined O_consthmix && !defined O_biharmonic && !defined O_isopycmix
# if defined O_bryan_lewis_horizontal
# else
# endif
#else
#endif
#if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
#endif
#if defined O_isopycmix
#endif
