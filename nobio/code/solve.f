! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers7_nobio/updates/solve.F
      subroutine solve (n)

!=======================================================================
!     solve for tracer distribution after diffusion

!     input:
!       n    = tracer number
!=======================================================================

      implicit none

      integer i, ii, ierr, j, jj, k, n

      logical done

      real afw, afe, afn, atc, ate, atn, ats, atnc, atsc, atw, b, dt
      real dtss, dfw, dfe, dfn, fa, fb, fc, fd, ff, fg, fh, tmp, x

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "solve.h"
      include "atm.h"
      include "cembm.h"
      include "csbc.h"
      include "grdvar.h"
      include "coord.h"
      include "levind.h"
      include "ice.h"

      real forc(imt,jmt)

      dtss = dts

!-----------------------------------------------------------------------
!     set the forcing for each tracer
!-----------------------------------------------------------------------

      if (n .eq. isat) then

!       temperature

        fa = dtss/(cpatm*rhoatm*sht)
        fb = dtss*vlocn/(cpatm*rhoatm*sht)
        fc = dtss*slice/(cpatm*rhoatm*sht) - fb
        fd = scatter*(1. + pass)
        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = fa*(solins(i,j)*sbc(i,j,iaca)*fd
     &                - dnswr(i,j)*scatter - outlwr(i,j)
     &                + uplwr(i,j) + upsens(i,j))
!           latent heat from total precipitation as water
            forc(i,j) = forc(i,j) + precip(i,j)*fb
!           correct for latent heat from snow
            forc(i,j) = forc(i,j) + fc*psno(i,j)
          enddo
        enddo

      elseif (n .eq. ishum) then

!       humidity

        fa = dtss/(rhoatm*shq)
        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = fa*evap(i,j)
          enddo
        enddo

      else

!       other tracers

        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = c0
          enddo
        enddo

      endif
      call embmbc (forc)

!-----------------------------------------------------------------------
!     calculate new coefficients if required
!-----------------------------------------------------------------------

      if (newcoef(lf,n)) call coef (n)

!-----------------------------------------------------------------------
!     shuffle in time
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          tmp = at(i,j,2,n)
          at(i,j,2,n) = at(i,j,lf,n) + forc(i,j)
          at(i,j,1,n) = tmp
        enddo
      enddo

!-----------------------------------------------------------------------
!     load rhs into the solver array
!-----------------------------------------------------------------------

      k = 0
      do j=2,jmtm1
        do i=2,imtm1
          k = k + 1
          bv(k) = at(i,j,2,n)
          xv(k) = at(i,j,1,n)
        enddo
      enddo

!-----------------------------------------------------------------------
!     solve for tracer
!-----------------------------------------------------------------------

      call mgrid (xv, ap(1,lf,n), an(1,lf,n), as(1,lf,n), ae(1,lf,n)
     &,           aw(1,lf,n), bv, 1, iimtm2, 1, jjmtm2, iimtm2, jjmtm2
     &,           itin(n), levelin, epsin(n), itout(n), levelout
     &,           epsout(n))

      newcoef(lf,n) = .false.
      if (epsout(n) .gt. epsin(n)) write(*,*)
     &  '==> Warning:  atmospheric solver not converging in ',
     &  itout(n),' iterations ( eps = ',epsout(n), ' > ',epsin(n),' )'

!-----------------------------------------------------------------------
!     copy new solution from left hand side
!-----------------------------------------------------------------------

      k = 0
      do j=2,jmtm1
        do i=2,imtm1
          k = k + 1
          at(i,j,2,n) = xv(k)
        enddo
      enddo

!-----------------------------------------------------------------------
!     set boundary conditions
!-----------------------------------------------------------------------

      call embmbc (at(1,1,2,n))

      return
      end

      subroutine coef (n)

!=======================================================================
!     compute matrix coefficients

!     input:
!       n    = tracer number
!=======================================================================

      implicit none

      integer i, ide, ii, ielm, iord, j, jdn, jj, jord, n, iwx, iwy

      real acej, acnj, adde, addn, adds, addw, cc, ce, cew, cmax, cn
      real cns, cs, cw, dcej, dcnj, fe, fn, fs, fw, ue, un, uw, ve
      real vn, vs, vw

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "solve.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "csbc.h"
      real diff_factor

      diff_factor = max(0., (1. + adiff*dtbar))

      ielm = 0
      iord = 0

      iwx = 0
      iwy = 0
      if (n .eq. isat) then
        iwx = iwxt
        iwy = iwyt
      elseif (n .eq. ishum) then
        iwx = iwxq
        iwy = iwyq
      elseif (n .eq. ico2) then
        iwx = iwxc
        iwy = iwyc
      endif

      if (iwx .gt. 0) call embmbc (sbc(1,1,iwx))
      if (iwy .gt. 0) call embmbc (sbc(1,1,iwy))

      cmax = 3.9e10

      do j=2,jmtm1
        jj = j - 1
        jdn = j
        do i=2,imtm1
          ii = i - 1
          ide = i

!-----------------------------------------------------------------------
!         set coefficients for implicit diffusion
!-----------------------------------------------------------------------

          cs = dn(i,j-1,n)
          cn = dn(i,jdn,n)
          cw = de(i-1,j,n)
          ce = de(ide,j,n)
          if (n .eq. isat) then
            cn = cn*diff_factor
            cs = cs*diff_factor
            ce = ce*diff_factor
            cw = cw*diff_factor
          endif

!-----------------------------------------------------------------------
!         closed north/south boundary conditions for diffusion
!-----------------------------------------------------------------------

          if (j .eq. 2) cs = c0
          if (j .eq. jmtm1) cn = c0

          cs =-dts*cs*dsgrd(j)
          cn =-dts*cn*dngrd(j)
          cw =-dts*cw*cstr(j)*cstr(j)*dwgrd(i)
          ce =-dts*ce*cstr(j)*cstr(j)*degrd(i)

          cc = 1.0 - cs - cn - cw - ce

!-----------------------------------------------------------------------
!         set coefficients for up-stream advection
!-----------------------------------------------------------------------

          vs = (sbc(i-1,j-1,iwy) + sbc(i,j-1,iwy))
          vn = (sbc(i-1,j,iwy) + sbc(i,j,iwy))
          uw = (sbc(i-1,j-1,iwx) + sbc(i-1,j,iwx))
          ue = (sbc(i,j-1,iwx) + sbc(i,j,iwx))

          fs = p5*(c1 + sign(c1,vs))
          fn = p5*(c1 + sign(c1,vn))
          fw = p5*(c1 + sign(c1,uw))
          fe = p5*(c1 + sign(c1,ue))

!-----------------------------------------------------------------------
!         closed north/south boundary conditions for advection
!-----------------------------------------------------------------------

          if (j .eq. 2) vs = c0
          if (j .eq. jmtm1) vn = c0

          cs = cs - dts*fs*vs*asgrd(j)
          cn = cn + dts*(c1-fn)*vn*angrd(j)
          cw = cw - dts*fw*uw*cstr(j)*azgrd(i)
          ce = ce + dts*(c1-fe)*ue*cstr(j)*azgrd(i)
          cc = cc + dts*(fn*vn*angrd(j)-(c1-fs)*vs*asgrd(j)
     &       + (fe*ue - (c1-fw)*uw)*cstr(j)*azgrd(i))

          iord = iord + 1

!-----------------------------------------------------------------------
!         load the coefficients for the multigrid solver
!-----------------------------------------------------------------------

          ap(iord,lf,n) = cc
          an(iord,lf,n) = -cn
          as(iord,lf,n) = -cs
          ae(iord,lf,n) = -ce
          aw(iord,lf,n) = -cw

        enddo
      enddo

      return
      end
