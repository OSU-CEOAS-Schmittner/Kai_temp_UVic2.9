gyre.f
size.h
param.h
pconst.h
stdunits.h
coord.h
cregin.h
diag.h
grdvar.h
hmixc.h
mw.h
scalar.h
isopyc.h
#if defined O_mom && defined O_gyre_components
# if defined O_isopycmix
# endif
# if defined O_consthmix && !defined O_biharmonic && !defined O_isopycmix
#  if defined O_bryan_lewis_horizontal
#  else
#  endif
# else
# endif
#endif
