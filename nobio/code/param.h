! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers7_nobio_linadv/updates/param.h
!======================= include file "param.h" ========================

!     main parameter file which sets ocean characteristics:

!     nvar   = number of prognostic variables
!     lseg   = maximum number of longitudinal stream function segments
!     nlatpr = maximum number of latitudes for matrix printouts
!              on diagnostic time steps
!     nhreg  = number of regions in the horizontal used for averaging
!              tracers.
!     nvreg  = number of regions in the vertical used for term balance
!              calculations. note "nvreg" is not used for tracer
!              averages
!     numreg = total number of regions ( = product of nhreg & nvreg)
!              used for term balance calculations

!     nvarbh = number of prognostic variables using biharmonic mixing

!     ncrows = number of calculated rows within the MW.
!              (the remaining rows are buffer rows).

      integer lseg, nlatpr, nhreg, nvreg, numreg, nvar, nvarbh
      integer imtm1, kmm1, imtp1, imtm2, jmtp1, jmtm1, jmtm2, jscan
      integer kmp1, kmp2, imtkm, nwds, nkflds, nslab, ntmin2, ncrows

      parameter (lseg=5, nlatpr=10)
      parameter (nhreg=3, nvreg=1, numreg=nhreg*nvreg)
      parameter (nvar=nt+2)

      parameter (imtm1=imt-1, kmm1=km-1)
      parameter (imtp1=imt+1, imtm2=imt-2
     &,          jmtp1=jmt+1, jmtm1=jmt-1, jmtm2=jmt-2
     &,          jscan=jmtm2
     &,          kmp1=km+1, kmp2=km+2
     &,          imtkm=imt*km, nwds=imt*jmt, nkflds=2
     &,          nslab=imt*nvar*km, ntmin2=nt+1/nt)

      parameter (ncrows = jmw - 3 + jmw/jmt)
