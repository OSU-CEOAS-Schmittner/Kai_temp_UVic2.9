! source file: /usr/local/models/UVic_ESCM/2.9/source/embm/solve.h
!======================== include file "solve.h" ========================

!     variables needed for solving atmospheric advection and diffusion

!     newcoef = logical flag for calculating new coefficients

      logical newcoef(2,nat)
      common /solver_l/ newcoef

      integer iimtm2, jjmtm2, nord, nelm
      parameter (iimtm2 = imtm2)
      parameter (jjmtm2 = jmtm2)
      parameter (nord = iimtm2*jjmtm2)

!     itin    = requested maximum iterations
!     itout   = actual iterations
!     newcoef = logical flag for calculating new coefficients
!     bv      = right hand side vector (b)
!     xv      = left hand side vector (x)
!     epsin   = requested maximum error
!     epsout  = actual error

      integer itin(nat), itout(nat)
      common /solver_i/ itin, itout

!     for mgrid routine storage is by compass coefficient
!     ap, an, as, ae, aw = centre, north, south, east, and west coef
!     ie. ap*xp = an*xn + as*xs + ae*xe + aw*xw + bp

      integer levelin, levelout
      common /solver_i/ levelin, levelout

      real(kind=8) bv, xv, epsin, epsout, an, as, ae, aw, ap
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solver_r/ an(nord,2,nat), as(nord,2,nat), ae(nord,2,nat)
      common /solver_r/ aw(nord,2,nat), ap(nord,2,nat)

!     grid terms for the atmospheric solver

      real dwgrd, degrd, azgrd, dsgrd, dngrd, asgrd, angrd
      common /solve_r/ dwgrd(2:imtm1), degrd(2:imtm1), azgrd(2:imtm1)
      common /solve_r/ dsgrd(2:jmtm1), dngrd(2:jmtm1), asgrd(2:jmtm1)
      common /solve_r/ angrd(2:jmtm1)
