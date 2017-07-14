! source file: /raid24/aho/UVic2.9/default_comb2/nobio/updates/hmixc.h
!======================= include file "hmixc.h" ========================

!                    horizontal mixing coefficients

!     visc_cnu = viscosity coeff for northern face of "u" cell
!     visc_ceu = viscosity coeff for eastern face of "u" cell
!     diff_cnt = diffusion coeff for northern face of "T" cell
!     diff_cet = diffusion coeff for eastern face of "T" cell

!     am     = constant lateral viscosity coeff for momentum
!     ah     = constant lateral diffusion coeff for tracers
!     am3    = viscosity coeff for metric term on "u" cell
!     am4    = another viscosity coeff for metric term on "u" cell
!     ambi   = constant lateral biharmonic viscosity coeff for momentum
!     ahbi   = constant lateral biharmonic diffusion coeff for tracers
!=======================================================================

      real am, ambi, am3, am4, ah, ahbi, visc_ceu, visc_cnu, amc_north
      real amc_south, Ahh(km), diff_cnt, diff_cet, ahc_north, ahc_south
      real strain, am_lambda, am_phi, smag_metric, diff_c_back
      real hl_depth, hl_back, hl_max, hl_u, hl_n, hl_e, hl_b
      real droz, rich_inv

      common /diffus_r/ am, ambi, am3(jmt), am4(jmt,2)
      common /diffus_r/ ah, ahbi
      common /diffus_r/ visc_ceu(imt,km,jmt)
      common /diffus_r/ visc_cnu(imt,km,jmt)
      common /diffus_r/ amc_north(imt,km,jmt)
      common /diffus_r/ amc_south(imt,km,jmt)
      common /diffus_r/ diff_cnt, diff_cet
      common /diffus_r/ ahc_north(jmt), ahc_south(jmt)
