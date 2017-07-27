mom_tavg.f
#if defined O_mom
# if defined O_save_npzd
# endif
# if defined O_save_npzd
# endif
# if defined O_units_time_years
#  if !defined O_save_time_relyear0
#  else
#  endif
# else
#  if !defined O_save_time_relyear0
#  else
#  endif
# endif
# if defined O_save_npzd
# endif
# if defined O_save_anisotropic_viscosity
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
# endif
# if defined O_stream_function
# endif
# if defined O_save_convection
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_save_npzd
#  if !defined O_npzd_caco3
#  endif
# endif
# if defined O_pipe_co2
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_gent_mcwilliams
# endif
# if defined O_save_kv
# endif
# if defined O_KGMdiag
# endif
# if defined O_save_npzd
# else
# endif
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif      
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_extra_diagnostics
#   if defined O_npzd_nitrogen      
#   endif      
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_o2      
#   endif      
#  endif
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_gent_mcwilliams
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_save_kv
# endif
# if defined O_KGMdiag
# endif
# if defined O_save_npzd
#  if defined O_save_carbon_totals
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif      
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_iron
#     if defined O_npzd_iron_diagnostics
#      if defined O_npzd_caco3
#      endif
#     endif
#  endif
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_save_convection
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
# endif
# if defined O_stream_function
# endif
# if defined O_pipe_co2
# endif
# if defined O_save_anisotropic_viscosity
# endif 
# if defined O_gent_mcwilliams
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_save_kv
# endif
# if defined O_KGMdiag
# endif
# if defined O_save_npzd
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  else      
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif      
#  endif
#  if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
#     if defined O_npzd_caco3
#     endif      
#    endif
#  endif
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_save_convection
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
# endif
# if defined O_stream_function
# endif
# if defined O_pipe_co2
# endif
# if defined O_save_anisotropic_viscosity
# endif
# if defined O_save_npzd
# endif
# if defined O_save_npzd
# endif
# if defined O_save_npzd
# endif
# if defined O_save_anisotropic_viscosity
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
# endif
# if defined O_stream_function
# endif
# if defined O_save_convection
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_save_npzd
#  if !defined O_npzd_caco3
#  endif
# endif
# if defined O_pipe_co2
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_units_temperature_Celsius
# else
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_gent_mcwilliams
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_save_kv
# endif     
# if defined O_KGMdiag
# endif
# if defined O_save_npzd
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif      
#  endif
#  if defined O_npzd_extra_diagnostics
#   if defined O_npzd_nitrogen      
#   endif      
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_o2      
#   endif      
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics 
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if !defined O_npzd_caco3      
#  endif      
# endif
# if defined O_save_npzd
# endif
#endif
