setvbc.f
size.h
param.h
pconst.h
stdunits.h
coord.h
csbc.h
grdvar.h
levind.h
scalar.h
mw.h
#if defined O_mom
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "mw.h"
#ifndef O_TMM
#else
#endif ! SPK O_TMM
# if defined O_gthflx
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
#    if defined O_npzd_nitrogen
#    endif
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
# endif
#ifndef O_TMM
#endif
#endif ! SPK no O_TMM
