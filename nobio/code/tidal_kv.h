! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers7_nobio/updates/tidal_kv.h
!===================== include file "tidal_kv.h" =======================

!     quantities for parameterization of tidal mixing

!     Simmons, H.L., S.R. Jayne, L.C. St. Laurent and A.J. Weaver, 2004:
!     Tidally driven mixing in a numerical model of the ocean general
!     circulation. Ocean Modelling, 6, 245-263.
!
!     Schmittner, A., and G. Egbert, 2013

      real zetar, ogamma, gravrho0r
      common /tdr/ zetar, ogamma, gravrho0r

!     edr = energy dissipation rate due to tides

      real edrm2,edrs2,edrk1,edro1
      common /tdr/ edrm2(imt,km,jmt),edrs2(imt,km,jmt)
      common /tdr/ edrk1(imt,km,jmt),edro1(imt,km,jmt)
