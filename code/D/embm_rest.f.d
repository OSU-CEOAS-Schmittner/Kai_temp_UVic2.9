embm_rest.f
size.h
param.h
pconst.h
stdunits.h
atm.h
cembm.h
coord.h
csbc.h
grdvar.h
ice.h
evp.h
tmngr.h
switch.h
iounit.h
npzd.h
#if defined O_embm
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_carbon_13 && defined O_carbon_13_coupled && !defined O_c13ccn_data && !defined O_c13ccn_data_transient
# endif
# if defined O_carbon && !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
# endif
# if defined O_carbon_13 && defined O_carbon_13_coupled && !defined O_c13ccn_data && !defined O_c13ccn_data_transient
# endif
# if defined O_carbon_14 && defined O_carbon_14_coupled && !defined O_c14ccn_data && !defined O_c14ccn_data_transient
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
# endif
# if defined O_ice_data_transient
# endif
# if defined O_sealev_data_transient
# endif
# if defined O_sealev || defined O_sealev_data
#  if defined O_sealev_data_transient
#  endif
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_carbon_13
# endif
# if defined O_ism || defined O_landice_data_transient
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_carbon_13
# endif
# if defined O_carbon_14
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
# endif
# if defined O_ism || defined O_landice_data
# endif
# if defined O_sealev || defined O_sealev_data
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_carbon_13
# endif
# if defined O_ism || defined O_landice_data
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_carbon_13
# endif
# if defined O_carbon_13
#  if !defined O_carbon_13_coupled
#  endif
# endif
# if defined O_carbon_14
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
# endif
# if defined O_ism || defined O_landice_data
# endif
# if defined O_sealev || defined O_sealev_data
# endif
# if defined O_carbon
#  if defined O_carbon_13
#  endif
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_nitrogen
#  if defined O_npzd_nitrogen_15
#  endif
# endif
# if defined O_npzd_iron
# endif
# if defined O_carbon_13
# endif
# if defined O_ism || defined O_landice_data
# endif
# if defined O_embm_awind || defined O_embm_adiff
# endif
# if defined O_embm_awind
# endif
# if defined O_ice
# endif
# if defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
# endif
# if defined O_ice_evp && defined O_ice
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
# endif
#endif
