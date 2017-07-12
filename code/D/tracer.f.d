tracer.f
size.h
param.h
pconst.h
stdunits.h
accel.h
coord.h
cregin.h
csbc.h
emode.h
grdvar.h
hmixc.h
levind.h
mw.h
scalar.h
switch.h
timeavgs.h
tmngr.h
vmixc.h
diaga.h
ice.h
npzd.h
atm.h
cembm.h
isopyc.h
fdift.h
ctavg.h
diag.h
iounit.h
#if defined O_mom
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "timeavgs.h"
#include "tmngr.h"
#include "vmixc.h"
# if defined O_save_convection || defined O_carbon_14
#include "diaga.h"
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
# if defined O_npzd || defined O_carbon_14
#include "npzd.h"
# endif
# if defined O_npzd
#  if defined O_npzd_o2
#  endif
#  if defined O_carbon_13 || defined O_npzd_caco3
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif     
#  if defined O_embm
#include "atm.h"
#   if defined O_carbon_13 || defined O_npzd_caco3
#include "cembm.h"
#   endif
#  endif
# endif
# if defined O_carbon_fnpzd
#include "calendar.h"
# endif
# if defined O_plume
# endif
# if defined O_npzd || defined O_carbon_14
#  ifdef O_TMM
#  endif 
# endif
# if defined O_carbon_fnpzd
# endif
# if defined O_isopycmix
#include "isopyc.h"
# endif
#include "fdift.h"
# if defined O_carbon_fnpzd
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
# endif
# ifndef O_TMM
#  if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic || defined O_isopycmix
#  else
#  endif
#  if defined O_consthmix && !defined O_bryan_lewis_horizontal && !defined O_biharmonic
#  endif
#  if defined O_plume
#  endif
# endif ! SPK not O_TMM
# if defined O_npzd
#  if defined O_ice
#   if defined O_ice_cpts
#   else
#   endif
#  else
#  endif
#  if defined O_embm
#    if defined O_npzd_caco3
#    else
#    endif
#  else
#    if defined O_npzd_caco3
#    else
#    endif
#  endif
#  if defined O_carbon            
#  endif            
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#    if defined O_npzd_caco3
#    endif            
#   endif
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
#   endif            
#  endif
#  if defined O_npzd_o2
#  endif
#  if defined O_carbon_13 || defined O_npzd_caco3
#   if defined O_carbon_co2_2d
#   else
#   endif
#   if defined O_carbon_13
#   endif
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#  if defined O_npzd_o2
#  endif
#  if defined O_carbon_13 || defined O_npzd_caco3
#  endif
#  if defined O_carbon_13
#  endif
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#  if defined O_sed
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_sed
#  endif
#  if defined O_time_averages && defined O_save_npzd
#   if !defined O_npzd_caco3
#   endif
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
#   if defined O_npzd_caco3
#   endif
#   if defined O_kk_ballast
#   endif              
#   if defined O_npzd_nitrogen && defined O_npzd_o2
#   endif
#   if defined O_npzd_iron
#   endif 
#  endif
# endif
# if defined O_carbon_fnpzd
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
# endif
# if defined O_carbon && defined O_carbon_14
#  if defined O_npzd
#  else
#  endif
# endif
# if defined O_sed && !defined O_sed_uncoupled
#  if defined O_carbon && defined O_npzd
#   if defined O_global_sums
#   endif
#  endif
#  if defined O_npzd_alk
#  endif
# endif
#ifndef O_TMM
# if defined O_consthmix
#  if !defined O_biharmonic || defined O_bryan_lewis_horizontal
#   if defined O_bryan_lewis_horizontal
#   else
#   endif
#   if defined O_isopycmix
#    if defined O_bryan_lewis_horizontal
#    else
#    endif
#   else
#   endif
#  else
#  endif
# else
#  if defined O_smagnlmix
#  endif
# endif
# if defined O_isopycmix
# endif
# if defined O_replacst
# else
# endif
# if defined O_source_term || defined O_npzd || defined O_carbon_14
#  if defined O_npzd || defined O_carbon_14
#  endif
#  if defined O_shortwave
#  endif
# endif
# if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
# endif
# if defined O_source_term || defined O_npzd || defined O_carbon_14
# endif
# if defined O_plume
# endif
# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
# endif
# if defined O_replacst
# endif
# if defined O_convect_brine
# else
#  if !defined O_implicitvmix || defined O_isopycmix
#   if defined O_fullconvect
#   else
#   endif
#  endif
# endif
# if defined O_save_convection
# endif
# if defined O_fourfil || defined O_firfil
# endif
# if !defined O_replacst
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
#  if !defined O_npzd_no_vflux
#   if defined O_kk_ballast
#   endif
#   if defined O_npzd_caco3
#   endif
#  endif
#  if defined O_npzd_nitrogen
#   if !defined O_npzd_no_vflux
#   endif
#   if defined O_npzd_nitrogen_15
#    if !defined O_npzd_no_vflux
#     if defined O_npzd_caco3
#     endif      
#    endif
#   endif
#  endif
#  if defined O_npzd_iron
#   if !defined O_npzd_no_vflux
#   endif
#  endif
#  if defined O_carbon_13
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_sed
#  if defined O_carbon
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_npzd_o2
#  endif
# endif
# if defined O_carbon && defined O_carbon_14
# endif
#endif ! SPK not O_TMM
#ifndef O_TMM
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
# if defined O_tracer_averages
#include "ctavg.h"
# endif
# if defined O_npzd_caco3
#include "npzd.h"
# endif      
#include "diag.h"
#include "diaga.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
# if defined O_meridional_tracer_budget
#include "ctmb.h"
# endif
# if defined O_time_step_monitor
# endif
# if defined O_isopycmix
#include "isopyc.h"
# endif
#include "fdift.h"
# if defined O_save_mixing_coeff
#  if !defined O_consthmix || defined O_biharmonic || defined O_isopycmix
#  else
#  endif
#  if defined O_isopycmix
#  endif
# endif
# if defined O_save_convection_full
# endif
# if defined O_time_step_monitor
# endif
# if defined O_tracer_averages
# endif
# if defined O_tracer_yz
#  if defined O_source_term || defined O_npzd || defined O_carbon_14
#  endif
# endif
# if defined O_meridional_tracer_budget
#  if defined O_source_term || defined O_npzd || defined O_carbon_14
#  endif
# endif
# if defined O_gyre_components
# endif
# if defined O_term_balances
# endif
# if defined O_xbts
# endif
# if defined O_mom_tbt
# endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "diaga.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "timeavgs.h"
# if defined O_save_convection_full
# endif
# if defined O_term_balances
# endif
# if defined O_xbts
# endif
# if defined O_mom_tbt
# endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "levind.h"
#include "mw.h"
#include "switch.h"
#include "vmixc.h"
#  if defined O_xbts || defined O_mom_tbt
#  else
#   if defined O_term_balances
#   endif
#  endif
#  if defined O_xbts || defined O_mom_tbt
#  else
#   if defined O_term_balances
#   endif
#  endif
# endif
# if defined O_mom && defined O_shortwave
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "cshort.h"
# endif
#endif ! SPK not O_TMM
#endif
