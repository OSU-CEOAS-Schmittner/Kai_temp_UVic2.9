! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/tbt.h
!====================== include file "tbt.h" ===========================

!     ntbtts = number of tracer term balance time steps

!     tbt   = time averagred term balances of tracers
!       (1) = total time rate of change for the tracer
!       (2) = change due to zonal nonlinear term (UT)x
!       (3) = change due to meridional nonlinear term (VT)y
!       (4) = change due to vertical nonlinear term (WT)z
!       (5) = change due to zonal diffusion: Ah*Txx
!       (6) = change due to meridional diffusion: Ah*Tyy
!       (7) = change due to vertical diffusion:  kappa_h*Tzz
!       (8) = change due to source term
!       (9) = change due to explicit convection
!      (10) = change due to filtering
!      (11) = average tracer within volume (tracer units)
