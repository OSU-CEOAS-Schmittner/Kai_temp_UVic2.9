embm_tsi.f
#if defined O_embm
# if defined O_units_time_years
#  if !defined O_save_time_relyear0
#  else
#  endif
# else
#  if !defined O_save_time_relyear0
#  else
#  endif
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_save_carbon_totals
# else
# endif
# if defined O_carbon_13
# endif
# if defined O_carbon_14
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_ice
#   if defined O_landice_data_transient || defined O_ism
#  endif
# endif
# if defined O_tai_ns
#  if defined O_units_temperature_Celsius
#  else
#  endif
#  if defined O_ice
#   if defined O_landice_data_transient || defined O_ism
#   endif
#  endif
# endif
# if defined O_tai_lo
#  if defined O_units_temperature_Celsius
#  else
#  endif
# endif
# if defined O_tai_rad
#  if defined O_tai_lo
#  endif
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_carbon_13
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_save_carbon_totals
# else
# endif
# if defined O_carbon_13
# endif
# if defined O_carbon_14
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_ice
#  if defined O_landice_data_transient || defined O_ism
#  endif
# endif
# if defined O_tai_ns
#  if defined O_units_temperature_Celsius
#  else
#  endif
#  if defined O_ice
#   if defined O_landice_data_transient || defined O_ism
#   endif
#  endif
# endif
# if defined O_tai_lo
#  if defined O_units_temperature_Celsius
#  else
#  endif
# endif
# if defined O_tai_rad
#  if defined O_tai_lo
#  endif
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_carbon_13
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
# endif
#endif
