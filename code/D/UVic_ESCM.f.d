UVic_ESCM.f
size.h
param.h
pconst.h
stdunits.h
coord.h
csbc.h
iounit.h
emode.h
levind.h
scalar.h
switch.h
tmngr.h
cembm.h
atm.h
mw.h
calendar.h
accel.h
cnep.h
cprnts.h
diag.h
fwa.h
hmixc.h
insolation.h
isopyc.h
mtlm.h
npzd.h
sed.h
stab.h
veg.h
vmixc.h
#ifndef O_TMM
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "iounit.h"
#include "emode.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "cembm.h"
#if defined O_embm
#include "atm.h"
#endif
#if defined O_mom
#include "mw.h"
#endif
#if defined O_mom
# if defined O_sed
# endif
#endif
#if defined O_embm
# if defined O_ism
# endif
# if defined O_mtlm
# endif
#endif
#if defined O_mtlm
# if defined O_mtlm_segday
# else
# endif
#endif
#if !defined O_embm && defined O_mom
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
#   if defined O_npzd_nitrogen
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
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
#    if defined O_npzd_caco3
#    endif            
#   if defined O_npzd_nitrogen
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#if defined O_embm
# if !defined O_mom
#  if defined O_mtlm
#  endif
# endif
#else
#endif
#if defined O_mtlm
#endif
#if defined O_mom
# if defined O_embm
# endif
# if defined O_ism
# endif
# if defined O_mtlm
# endif
# if defined O_sed
# endif
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#endif ! O_TMM
#ifndef O_TMM
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "cembm.h"
#include "calendar.h"
#include "switch.h"
#include "tmngr.h"
#if defined O_sealev_salinity && defined O_mom
#include "mw.h"
#include "iounit.h"
#include "levind.h"
#include "grdvar.h"
#endif
# if defined O_even_fluxes
#endif
#if defined O_embm && !defined O_mom && !defined O_replacst
#endif
#if defined O_redi_diffusion || defined O_gent_mcwilliams
# if !defined O_isopycmix
# endif
#endif
#if defined O_mom
#endif
#if defined O_embm
#endif
#if defined O_ism
#endif
#if defined O_mtlm
#endif
#if defined O_sed
#endif
#if defined O_mom && defined O_embm && defined O_restorst
#endif
#if defined O_mom && defined O_embm && defined O_replacst
#endif
#if defined O_restorst && defined O_replacst
#endif
#if defined O_co2ccn_data_transient && !defined O_co2ccn_data
#endif
#if defined O_co2emit_data_transient && !defined O_co2emit_data
#endif
#if defined O_co2emit_track_sat_transient && !defined O_co2emit_track_sat
#endif
#if defined O_co2emit_track_co2_transient && !defined O_co2emit_track_co2
#endif
#if defined O_agric_data_transient && !defined O_agric_data
#endif
#if defined O_crop_data_transient && !defined O_crop_data
#endif
#if defined O_landice_data_transient && !defined O_landice_data
#endif
#if defined O_solar_data_transient && !defined O_solar_data
#endif
#if defined O_volcano_data_transient && !defined O_volcano_data
#endif
#if defined O_sulphate_data_transient && !defined O_sulphate_data
#endif
#if defined O_aggfor_data_transient && !defined O_aggfor_data
#endif
#if defined O_cfcs_data_transient && !defined O_cfcs_data
#endif
#if defined O_c14ccn_data_transient && !defined O_c14ccn_data
#endif
#if defined O_co2ccn_data|| defined O_co2emit_data || defined O_co2emit_track_co2 || defined O_co2emit_track_sat
# if defined O_co2ccn_data
# endif
# if defined O_co2emit_data
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat
# endif
# if defined O_co2ccn_user
# endif
#endif
#endif ! O_TMM
#ifndef O_TMM
#if defined O_mom
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "iounit.h"
#include "mw.h"
#include "tmngr.h"
#endif
#endif ! O_TMM
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#if defined O_embm_awind && defined O_embm
#endif
#if defined O_embm
# if defined O_carbon_co2_2d
# endif
#endif
#if defined O_shortwave
#endif
#if defined O_ice_evp
#endif
#if defined O_carbon
# if defined O_carbon_13
# endif
# if defined O_carbon_14
# endif
#endif
#if defined O_npzd_alk
#endif
#if defined O_npzd_o2
#endif
#if defined O_npzd
# if !defined O_npzd_no_vflux
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_iron
#  if !defined O_npzd_no_vflux
#  endif
# endif
# if defined O_npzd_nitrogen
#  if !defined O_npzd_no_vflux
#  endif
#  if defined O_npzd_nitrogen_15
#   if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
# endif
# if defined O_carbon_13
#  if !defined O_npzd_no_vflux
#    if defined O_npzd_caco3
#    endif 
#   if defined O_npzd_nitrogen
#   endif
#  endif
# endif
#endif
#if defined O_cfcs_data || defined O_cfcs_data_transient
#endif
#if defined O_mtlm
#endif
#if defined O_mtlm && defined O_carbon
#endif
#if defined O_mtlm_carbon_13
#endif
#if defined O_mtlm_carbon_14
#endif
#if defined O_sed
# if defined O_carbon
# endif
# if defined O_npzd_alk
# endif
#endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "atm.h"
#if defined O_mom
#include "mw.h"
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd
#  if defined O_npzd_caco3
#  endif     
#  if defined O_kk_ballast
#  endif       
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
#   endif      
#   if defined O_npzd_nitrogen
#   endif
#  endif     
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_npzd
#  if defined O_carbon
#   if defined O_carbon_13
#   endif
#  endif
#  if defined O_npzd_alk
#  endif      
#  if defined O_npzd_o2
#  endif
#  if defined O_npzd_iron
#  endif
#  if defined O_kk_ballast
#  endif      
#  if defined O_npzd_caco3
#  endif      
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
#   endif      
#   if defined O_npzd_nitrogen
#   endif
#  endif
# endif
#endif
#if defined O_embm
# if defined O_carbon && defined O_carbon_co2_2d
# endif
#endif
#if !defined O_mom || !defined O_time_step_monitor
#endif
#if !defined O_mom || !defined O_isopycmix
#endif
#if !defined O_mom || !defined O_isopycmix || !defined O_gent_mcwilliams
#endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "calendar.h"
#include "cembm.h"
#include "coord.h"
#include "cnep.h"
#include "csbc.h"
#include "cprnts.h"
#include "diag.h"
#include "emode.h"
#include "fwa.h"
#include "hmixc.h"
#include "insolation.h"
#include "iounit.h"
#if defined O_mom && defined O_isopycmix
#include "isopyc.h"
#endif
#include "mtlm.h"
#include "npzd.h"
#include "scalar.h"
#include "sed.h"
#include "stab.h"
#include "switch.h"
#include "tmngr.h"
#include "veg.h"
#include "vmixc.h"
# if defined O_npzd_iron
# endif
#if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
#else
#endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_implicitvmix
# else
# endif
# if defined O_embm_snow_transient
# endif
# if defined O_embm
# else
# endif
# if defined O_mom
#  if defined O_isopycmix
#  endif
# else
# endif
# if defined O_carbon_14
# endif
# if defined O_carbon_13
# endif
# if defined O_calendar_360_day
# elif defined O_calendar_gregorian
# else
# endif
