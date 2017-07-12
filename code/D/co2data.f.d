co2data.f
size.h
param.h
pconst.h
stdunits.h
calendar.h
cembm.h
switch.h
tmngr.h
atm.h
#if defined O_co2emit_data || defined O_co2emit_data_transient
# if defined O_co2emit_data_transient
# else
# endif
# if defined O_co2emit_data_fuel && !defined O_co2emit_data_land
# elif defined O_co2emit_data_land && !defined O_co2emit_data_fuel
# else
# endif
#endif
#if defined O_co2ccn_data || defined O_co2ccn_data_transient || defined O_co2emit_track_co2 || defined O_co2emit_track_co2_transient
# if defined O_carbon_co2_2d
# endif
# if defined O_co2emit_track_co2 && !defined O_co2emit_track_sat
# endif
# if defined O_co2emit_track_co2
# else
# endif
# if defined O_co2ccn_data_transient || defined O_co2emit_track_co2_transient
# else
# endif
# if defined O_co2emit_track_co2
#  if defined O_carbon_co2_2d
#  endif
#  if !defined O_co2emit_track_sat
#  endif
# endif
# if defined O_co2emit_track_sat
#  if defined O_co2emit_track_co2
#  else
#  endif
# endif
# if defined O_co2ccn_data || defined O_co2ccn_data_transient
#  if defined O_carbon_co2_2d
#  endif
# else
#  if defined O_carbon_co2_2d
#  endif
# endif
#endif
#if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient || defined O_embm_vcs
# if defined O_landice_data
#  if defined O_ice_cpts
#  endif
# endif
# if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient
# if defined O_co2emit_track_sat_transient
# else
# endif
# endif
# if defined O_landice_data
# endif
# if defined O_sealev || defined O_sealev_data
# endif
# if defined O_co2emit_track_co2
# endif
# if defined O_co2emit_track_sat
# endif
#endif
