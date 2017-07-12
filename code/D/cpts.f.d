cpts.f
#if defined O_ice && defined O_ice_cpts && defined O_embm
# if defined O_ice_cpts_warnings
# endif
# if defined O_ice_evp
# endif
# if defined O_ice_cpts_simple_growth
# else
# endif
# if defined O_ice_cpts_roth_press
# endif
# if defined O_ice_cpts_roth_press
# endif
# if defined O_ice_cpts_roth_press
# endif
# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
# endif
# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
# endif
# if defined O_ice_cpts_roth_press
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# else
# endif
# if defined O_sealev || defined O_sealev_data
# endif
# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
# endif
# if defined O_ice_cpts_warnings
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# else
# endif
# if defined O_convect_brine
# else
# endif
# if defined O_convect_brine
# else
# endif
# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
# endif
# if defined O_ice_cpts_warnings
#  if !defined O_convect_brine
#  endif
# endif
#endif
