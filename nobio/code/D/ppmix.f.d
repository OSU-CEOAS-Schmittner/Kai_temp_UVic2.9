ppmix.f
#if defined O_mom && defined O_ppvmix
# if defined O_implicitvmix
# else
# endif
# if defined O_ppvmix && !defined O_implicitvmix
#  if defined O_isopycmix
#  endif
# endif
# if !defined O_implicitvmix
# else
# endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_implicitvmix
# else
# endif
#if defined O_bryan_lewis_vertical
#endif
# if defined O_matrix_sections
# endif
# if defined O_matrix_sections
# endif
#endif
