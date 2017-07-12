! source file: /raid24/aschmitt/UVic2.9/karin/mwc15_npzd_fe_n15_c13_alk_caco3/updates/size.h
!======================= include file "size.h" =========================

!-----------------------------------------------------------------------
!     USER INPUT:
!-----------------------------------------------------------------------

!     imt    = number of grid points in the longitudinal direction
!              (calculated points are from 2 through imt-1. End points
!               are boundaries)
!     jmt    = number of grid points (latitude rows) in the latitudinal
!              direction (calculated points are from 2 through jmt-1.
!              End points are boundaries)
!     km     = number of grid points in the vertical direction
!              (calculated points are from 1 through km)
!     nt     = number of tracers (temperature, salinity, ...)
!     nsrc   = number of tracer with sources
!     kpzd   = depth for limited npzd model
!     jmz    = size for "unrotated" zonal averages
!     jmzm1  = jmz minus one
!     mnisle = maximum number of islands (unconnected land masses)
!     maxipp = maximum number of all island perimeter points
!-----------------------------------------------------------------------

      integer imt, jmt, km, nt, nsrc, kpzd, nat, jmz, jmzm1, mnisle
      integer maxipp, jmw, jsmw, jemw
      parameter (imt=  102, jmt=  102, km= 19)
      parameter (nt=2 ! temp, salt
     $             +1 ! dic
     $             +1 ! dic13
     $             +1 ! alk
     $             +1 !
     $             +4 ! po4, phyt, zoop, detr
     $             +2 ! cocc, caco3
     $             +4 ! no3, diaz, dop, don
     $             +6 ! din15, phytn15, zoopn15, detrn15, diazn15, don15
     $             +1 ! coccn15
     $             +3 ! phytc13, zoopc13, detrc13
     $             +2 ! coccc13, caco3c13
     $             +2 ! diazc13, doc13
     $               +2 ! dfe, detrfe
     $               )
      parameter (nsrc=0
     $               +1 ! dic
     $               +1 ! dic13
     $               +1 ! alk
     $               +1 ! o2
     $               +4 ! po4, phyt, zoop, detr
     $               +2 ! cocc, caco3
     $               +4 ! no3, diaz, dop, don
     $               +6 ! din15, phytn15, zoopn15, detrn15, diazn15, don15
     $               +1 ! coccn15
     $               +3 ! phytc13, zoopc13, detrc13
     $               +2 ! coccc13, caco3c13
     $               +2 ! diazc13, doc13
     $               +2 ! dfe, detrfe
     $                 )
      parameter (kpzd=km)

      parameter (nat=2
     $, jmz=jmt, jmzm1=jmz-1)
      parameter (mnisle=50, maxipp=5000)

      parameter (jmw=jmt)

!-----------------------------------------------------------------------
!     set first and last calculated row within the MW. other rows
!     are used as buffers
!-----------------------------------------------------------------------

!     jsmw   = 1st calculated row within the MW
!     jemw   = last calculated row within the MW
      parameter (jsmw=2, jemw=jmw-1)
! Moses-Triffid land model

! POINTS = Maximum number of points in grid.
! STEPSM = Maximum number of timesteps in a day.
! klmax = maximum ocean depth levels over which the land model can exist

      integer POINTS, STEPSM, klmax
      parameter (POINTS=14300, STEPSM=24, klmax=0)

! NNVG  = Number of non-vegetation surface types.
! NPFT  = Number of plant functional types.
! NTYPE = Number of surface types.
! SOIL  = Index of the surface type 'Soil'
! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Soil

      integer NNVG, NPFT, NTYPE, SOIL
      parameter (NNVG=4, NPFT=5, NTYPE=6, SOIL=6)
