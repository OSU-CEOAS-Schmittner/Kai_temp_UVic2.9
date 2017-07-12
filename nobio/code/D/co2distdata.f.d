co2distdata.f
#if defined O_carbon_co2_2d
# if defined O_co2emit_data || defined O_co2emit_data_transient
# if defined O_co2emit_data_fuel && !defined O_co2emit_data_land
# elif defined O_co2emit_data_land && !defined O_co2emit_data_fuel
# else
# endif
#  if defined O_co2emit_data_transient
#  if defined O_co2emit_data_transient_repyr
#  endif
#  else
#  endif
# endif
#endif
