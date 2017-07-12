npzd.h
#if defined O_carbon
# if defined O_carbon_13
#  if defined O_npzd_nitrogen
#  endif		 
#  if defined O_npzd_caco3
#  endif		 
# endif
#endif		 
#if defined O_npzd_nitrogen
# if defined O_npzd_nitrogen_15
#  if defined O_npzd_caco3
#  endif		 
# endif
#endif
#if defined O_npzd_caco3
#endif
#if defined O_kk_ballast
#endif
#if defined O_npzd_iron
#endif
#if defined O_carbon
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
#endif
#if defined O_npzd_nitrogen
# if defined O_npzd_nitrogen_15
#  if defined O_npzd_caco3
#  endif
# endif
#endif
#if defined O_npzd_caco3
#endif
#if defined O_kk_ballast
#endif
#if defined O_npzd_iron
#endif
#if defined O_npzd_nitrogen_15
#endif
#if defined O_carbon_13
#endif
#if defined O_carbon_14
#endif
#if defined O_npzd
# if defined O_npzd_iron
# endif
# if defined O_npzd_iron
# endif
# if defined O_save_npzd
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_extra_diagnostics
#  endif
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics
#    if defined O_npzd_caco3
#    endif     
#   endif
#  endif
#  if defined O_npzd_nitrogen
#  endif
# endif
#endif
