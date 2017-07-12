npzd_src.f
size.h
npzd.h
calendar.h
coord.h
stdunits.h
scalar.h
param.h
pconst.h
grdvar.h
mw.h
#if defined O_mom && defined O_npzd 
#include "size.h"
#include "npzd.h"
#include "calendar.h"
#include "coord.h"
#include "stdunits.h"
#include "scalar.h"
# if defined O_npzd_iron
# endif
# if defined O_npzd_nitrogen
#  if !defined O_npzd_o2
#  endif      
#  if defined O_npzd_caco3
#  else      
#  endif      
# else
#  if defined O_npzd_caco3
#  else      
#  endif      
# endif
# if defined O_npzd_nitrogen_15
#  if !defined O_npzd_nitrogen
#  endif      
# endif
# if defined O_npzd_iron
#  if !defined O_npzd_o2
#  endif      
#  if defined O_npzd_caco3
#  endif      
# endif
# if defined O_npzd_caco3
#  if !defined O_carbon
#  endif      
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_caco3
# endif            
# if defined O_npzd_caco3
# endif 
# if defined O_npzd_nitrogen
#  if !defined O_npzd_caco3
#  else
#  endif
# else
#  if !defined O_npzd_caco3
#  else
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_carbon || defined O_npzd_alk
# endif
# if defined O_carbon
#  if defined O_carbon_13
#   if defined O_npzd_nitrogen
#   endif
#   if defined O_npzd_caco3
#   endif
#  endif
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif
#  endif
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_o2
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
# endif
# if defined O_carbon_13      
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_sed
# endif
#include "size.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "npzd.h"
# if defined O_npzd_extra_diagnostics
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd_iron
#  if defined O_npzd_iron_diagnostics
#  endif      
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
# if defined O_npzd_caco3
# endif
# endif
# if defined O_npzd_nitrogen
# endif   
# if defined O_sed
# endif
# if defined O_save_npzd
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
#  if defined O_carbon_13
#   if defined O_npzd_caco3         
#   endif
#  endif
#  if defined O_npzd_caco3         
#  endif         
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_nitrogen
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_o2
# else
# endif
# if defined O_npzd_nitrogen         
# endif
# if defined O_npzd_caco3
# endif        
# if defined O_kk_ballast
# endif
# if defined O_npzd_nitrogen
# endif
# if defined O_carbon
# endif
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif
#  endif
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_iron
#  if defined O_npzd_iron_diagnostics
#   if defined O_npzd_caco3
#   endif         
#  endif
# endif         
# if defined O_npzd_o2
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_caco3
# endif
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
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
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_npzd_caco3
#  if defined O_sed
#  endif
# endif
# if defined O_sed
#  if defined O_kk_ballast                   
#  else
#  endif
# endif
# if defined O_kk_ballast
# else               
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  else
#  endif
#  if defined O_kk_ballast
#  else               
#  endif
#  if defined O_save_npzd
#  endif
#  if defined O_npzd_nitrogen_15
#  endif              
# endif
# if defined O_npzd_iron
#  if defined O_npzd_iron_diagnostics
#  endif
# endif
# if defined O_carbon
# endif
# if defined O_carbon_13
# endif
# if defined O_save_npzd
#  if defined O_kk_ballast
#  else
#  endif               
# endif
# if defined O_carbon         
# endif         
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif         
#  endif
# endif
# if defined O_carbon_13
#  if defined O_npzd_nitrogen
#  endif         
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
# endif
# if !defined O_carbon
# else
#  if !defined O_npzd_caco3         
#  endif         
#  if defined O_carbon_13
#   if defined O_npzd_caco3          
#   endif
#   if !defined O_npzd_caco3          
#   endif
#  endif
# endif       
# if defined O_npzd_alk
#  if !defined O_npzd_caco3
#  endif
# endif
# if defined O_kk_ballast
# endif         
# if defined O_npzd_iron
# endif
# if defined O_npzd_caco3
# endif
# if defined O_save_npzd
# endif
# if defined O_npzd_o2
#  if defined O_npzd_nitrogen
#   if defined O_npzd_nitrogen_15
#   else
#   endif
#   if defined O_npzd_nitrogen_15
#   endif
#   if defined O_npzd_nitrogen_15
#   endif
#   if defined O_npzd_alk
#   endif
#   if defined O_save_npzd
#   endif
#  endif
# endif
# if defined O_carbon
#  if defined O_npzd_caco3
#  else
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
#   else              
#   endif
#  endif
# endif
# if defined O_npzd_alk
#   if defined O_npzd_caco3
#   else
#   endif
# endif
# if defined O_sed
# endif
# if defined O_carbon
#  if defined O_npzd_caco3
#  else
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
#   else  
#   endif
#  endif
# endif
# if defined O_npzd_alk
#  if defined O_npzd_caco3
#  else
#  endif
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_nitrogen
# endif      
# if defined O_carbon
# endif      
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_iron
#  if defined O_npzd_iron_diagnostics
#   if defined O_npzd_caco3
#   endif  
#  endif
# endif       
# if defined O_npzd_o2
# endif
#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "npzd.h"
# if defined O_npzd_iron
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_nitrogen
#  endif
# endif
# if defined O_carbon      
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif      
#  endif
# endif
# if defined O_carbon_13
#  if defined O_npzd_nitrogen
#  endif
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif      
#  endif
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif      
#  if defined O_npzd_nitrogen
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_carbon
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif      
#  endif
# endif
# if defined O_carbon_13
#   if defined O_npzd_caco3
#   endif 
#   if defined O_npzd_nitrogen
#   endif 
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_nitrogen
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_iron
# else
# endif
# if defined O_npzd_iron
# else
# endif
# if defined O_npzd_nitrogen
# endif      
# if defined O_npzd_caco3
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_iron
#  else
#  endif
#  if defined O_npzd_iron
#  else
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_caco3
#  if defined O_npzd_iron
#  else
#  endif
#  if defined O_npzd_iron
#  else
#  endif
#  if defined O_npzd_nitrogen
#  endif
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_nitrogen
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
# endif
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#    if defined O_npzd_caco3
#    endif
#  endif
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics
#    if defined O_npzd_caco3
#    endif      
#   endif
#  endif
#  if defined O_npzd_extra_diagnostics
#  endif
# endif
# if defined O_npzd_iron
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_nitrogen
#  endif
# endif 
# if defined O_npzd_nitrogen
#  if defined O_npzd_caco3
#  endif        
# else
#  if defined O_npzd_caco3
#  endif
# endif        
# if defined O_npzd_nitrogen
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  endif
# else
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  endif
#  if defined O_npzd_caco3
#  endif        
# endif
# if defined O_npzd_caco3
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_nitrogen        
# endif        
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_iron
#  if defined O_npzd_nitrogen        
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif        
#  if defined O_kk_ballast
#  endif          
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_caco3
#   if defined O_npzd_nitrogen_15
#   endif
#  endif
#   if defined O_npzd_nitrogen_15
#   endif
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_nitrogen_15
#  endif
# else
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_kk_ballast
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_npzd_caco3
#  if defined O_npzd_nitrogen_15
#  endif        
#  if defined O_npzd_nitrogen_15
#  endif
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif        
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif        
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif  
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif
# if defined O_npzd_caco3
# endif
# if defined O_kk_ballast
# endif  
# if defined O_npzd_nitrogen
# endif
# if defined O_npzd_nitrogen_15
#  if defined O_npzd_caco3
#  endif        
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_nitrogen
#  endif        
# endif
# if defined O_npzd_caco3
#  if defined O_kk_ballast
#  endif
# else
# endif        
# if defined O_npzd_nitrogen
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_caco3
#   endif
#  else        
#   if defined O_npzd_caco3
#   endif
#  endif
#  if defined O_carbon
#   if defined O_npzd_caco3
#   endif
#  endif
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif
# else
#  if defined O_npzd_caco3
#  endif
#  if defined O_kk_ballast
#  else        
#   if defined O_npzd_caco3
#   endif
#  endif        
#  if defined O_carbon
#   if defined O_npzd_caco3
#   endif
#  endif
# endif
# if defined O_npzd_caco3
# endif        
# if defined O_npzd_iron
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_caco3
#   endif        
#  else
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_caco3
#   endif
#   if defined O_kk_ballast
#   endif        
#  endif
# endif
# if defined O_npzd_nitrogen_15
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_caco3
#  endif
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_caco3
#  endif        
# endif
# if defined O_carbon_13
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
#   endif        
#   if defined O_npzd_caco3
#   endif       
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_caco3
#   endif        
#   if defined O_npzd_caco3
#   endif        
#  else ! no nitrogen
#   if defined O_npzd_caco3
#   endif        
#   if defined O_npzd_caco3
#   endif
#   if defined O_npzd_caco3
#   endif        
#   if defined O_kk_ballast
#   else
#   endif        
#   if defined O_npzd_caco3
#   endif        
#  endif
# endif
# if defined O_npzd_nitrogen_15
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif
# endif
# if defined O_npzd_nitrogen
# endif
# if defined O_npzd_iron
# endif
# if defined O_save_npzd
#  if defined O_npzd_nitrogen
#   if defined O_npzd_caco3
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
#  if defined O_npzd_extra_diagnostics
#   if defined O_npzd_nitrogen
#   endif
#  endif
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif        
#  endif
# endif
# if defined O_npzd_caco3
#  if defined O_kk_ballast
#  endif
# endif        
# if defined O_npzd_iron
# endif
# if defined O_carbon_13
#  if defined O_npzd_caco3
#  endif        
#  if defined O_npzd_nitrogen
#  endif        
# endif
# if defined O_carbon
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#   if defined O_npzd_caco3
#   endif
#  endif
# endif
# if defined O_npzd_caco3
#  if defined O_kk_ballast
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_carbon_13
#  if defined O_npzd_nitrogen
#  endif      
#  if defined O_npzd_caco3
#  endif
# endif
#endif
