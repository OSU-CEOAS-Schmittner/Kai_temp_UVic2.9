gasbc.f
size.h
npzd.h
param.h
pconst.h
stdunits.h
coord.h
csbc.h
mw.h
ice.h
switch.h
tmngr.h
cembm.h
atm.h
insolation.h
calendar.h
grdvar.h
levind.h
solve.h
mtlm.h
#if defined O_embm
#include "size.h"
#include "npzd.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#  if defined O_mom
#include "mw.h"
#  endif
#  if defined O_ice
#   if defined O_ice_cpts
#include "cpts.h"
#   endif
#include "ice.h"
#  endif
#include "switch.h"
#include "tmngr.h"
#include "cembm.h"
#include "atm.h"
#include "insolation.h"
#include "calendar.h"
#include "grdvar.h"
#include "levind.h"
#include "solve.h"
#  if defined O_mtlm
#include "mtlm.h"
#  endif
#  if defined O_save_carbon_carbonate_chem
#include "diaga.h"
#  endif
# if !defined O_embm_annual
# endif
#ifndef O_TMM
#else
#endif
# if defined O_carbon
#  if defined O_carbon_13
#   if defined O_carbon_13_coupled
#   else
#   endif
#  endif
# endif
# if defined O_mom
#  if !defined O_constant_flux_reference || defined O_cfcs_data || defined O_cfcs_data_transient
#  endif
#  if !defined O_constant_flux_reference
#   if defined O_npzd_iron
#   endif
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_iron
#    endif
#    if defined O_kk_ballast
#    endif
#    if defined O_npzd_caco3
#    endif
#   endif
#   if !defined O_npzd_no_vflux
#   endif
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif      
#   endif
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if defined O_cfcs_data || defined O_cfcs_data_transient
#  endif
# endif
# if defined O_plume
# endif
# if defined O_convect_brine
# endif
# if defined O_mtlm
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
#  if defined O_npzd_iron
#   if !defined O_npzd_no_vflux
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
#  if defined O_carbon_13
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif
#    if defined O_npzd_nitrogen
#    endif
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_solar_data || defined O_solar_data_transient
# endif
# if !defined O_embm_annual
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
# endif
# if defined O_co2emit_data_transient
#  if defined O_carbon_co2_2d
#  endif
# endif
# if defined O_co2ccn_data || defined O_co2ccn_data_transient || defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
# endif
# if defined O_carbon_14
#  if defined O_c14ccn_data || defined O_c14ccn_data_transient
#   if defined O_c14ccn_data
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
# endif
# if defined O_embm_awind
# endif
# if defined O_embm && defined O_sealev_data_transient
#  if defined O_sealev_data_transient &&  defined O_sealev_salinity
#  else
#  endif
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_mom
#  if defined O_carbon || defined O_npzd_o2 || defined O_cfcs_data || defined O_cfcs_data_transient
#   if defined O_ice
#    if defined O_ice_cpts
#    else
#    endif
#   else
#   endif
#  endif
#  if defined O_carbon
#   if defined O_npzd_alk
#   else
#   endif
#   if defined O_carbon_co2_2d
#   else
#   endif
#   if defined O_save_carbon_carbonate_chem
#   endif
#   if defined O_carbon_co2_2d
#   endif
#   if defined O_carbon_13
#   endif
#   if defined O_carbon_14
#    if defined O_c14ccn_data
#    endif
#   endif
#  endif
#  if defined O_npzd_o2
#  endif
#  if defined O_cfcs_data || defined O_cfcs_data_transient
#  endif
# endif
# if defined O_carbon && defined O_mtlm
#  if defined O_carbon_co2_2d
#  else
#   if defined O_carbon_13
#    if defined O_mtlm_carbon_13
#    else
#    endif
#   endif
#  endif
#  if defined O_carbon_14
#   if defined O_mtlm_carbon_14
#   else
#   endif
#  endif
# endif
# if defined O_carbon
#  if defined O_mtlm && defined O_global_sums
#  endif
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_co2_2d
#   if !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
#   endif
#  else
#   if !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
#   endif
#   if defined O_carbon_13
#    if defined O_carbon_13_coupled
#    endif
#   endif
#  endif
#  if defined O_global_sums
#  endif
#  if defined O_carbon_14
#   if defined O_carbon_14_coupled
#   endif
#   if defined O_c14ccn_data
#   endif
#  endif
# endif
# if defined O_npzd_o2
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_crop_data_transient || defined O_pasture_data_transient || defined O_agric_data_transient
# endif
# if defined O_time_averages
# endif
# if defined O_time_step_monitor
# endif
#endif
