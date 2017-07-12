timeavgs.f
size.h
stdunits.h
param.h
pconst.h
levind.h
cregin.h
timeavgs.h
npzd.h
coord.h
docnam.h
csbc.h
emode.h
grdvar.h
iounit.h
scalar.h
switch.h
tmngr.h
calendar.h
cembm.h
atm.h
hmixc.h
#if defined O_mom
# if defined O_time_averages
#  if defined O_embm
#  endif
#  if defined O_save_anisotropic_viscosity
#  endif
#  if defined O_embm && !defined O_save_fluxes_with_virtual
#  else
#  endif
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  endif
#  if defined O_stream_function
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#    if defined O_npzd_caco3
#    endif            
#   endif
#   if defined O_kk_ballast
#   endif
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
#     if defined O_npzd_caco3
#     endif            
#    endif
#   endif
#   if defined O_npzd_extra_diagnostics
#   endif
#   if defined O_kk_ballast
#   endif
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_nitrogen
#   endif
#  if !defined O_npzd_caco3           
#  endif
#  endif
#  if defined O_save_convection
#  endif
# if defined O_pipe_co2
# endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if !defined O_save_time_relyear0
#  endif
#  if defined O_save_time_endper
#  elif defined O_save_time_startper
#  else
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  else
#   if defined O_calendar_360_day
#   elif defined O_calendar_gregorian
#   else
#   endif
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_save_carbon_totals
#   endif
#    if !defined O_npzd_caco3 
#    else
#    endif        
#   if defined O_kk_ballast
#   endif
#   if defined O_npzd_nitrogen
#    if defined O_npzd_caco3
#    endif        
#   endif
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
#     if defined O_npzd_caco3
#     endif
#    endif
#   endif
#   if defined O_npzd_extra_diagnostics
#   endif 
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  endif
#  if defined O_stream_function
#  endif
# if defined O_pipe_co2
# endif
#  if defined O_save_anisotropic_viscosity
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#    if defined O_npzd_caco3
#    endif      
#   endif
#   if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
#     if defined O_npzd_caco3
#     endif      
#    endif
#   endif
#   if defined O_kk_ballast
#   endif
#   if defined O_npzd_caco3
#   else
#   endif
#   if defined O_npzd_extra_diagnostics
#   endif
#  endif
#  if defined O_save_convection
#  endif
# if defined O_pipe_co2
# endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#    if defined O_npzd_caco3
#    endif     
#   endif
#   if defined O_npzd_iron
#    if defined O_npzd_iron_diagnostics
#     if defined O_npzd_caco3
#     endif
#    endif
#   endif
#   if defined O_npzd_caco3
#   else
#    if defined O_kk_ballast
#    endif
#   endif
#   if defined O_npzd_extra_diagnostics
#   endif
#  endif
#  if defined O_save_convection
#  endif
# if defined O_pipe_co2
# endif
#  if defined O_save_carbon_carbonate_chem
#  endif
# endif
#endif
