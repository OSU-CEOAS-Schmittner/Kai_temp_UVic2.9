! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers7_nobio/updates/hmixc.F
      subroutine hmixc (joff, js, je, is, ie)

!=======================================================================
!     set horizontal mixing coeffs on north and east face of "t" and
!     "u" cells.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer js ,je ,jrowstart ,joff ,jrowend ,jrow ,jm1 ,jp1 ,k ,ie
      integer is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "atm.h"
      include "levind.h"
      include "coord.h"

      integer i, istrt, iend, j, jstrt, jend
      real aeddy, beddy, beta, delx, WBC, gridlen, V0, Vz, efD
      real px, N, bmunk

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     set all horizontal mixing coefficients
!-----------------------------------------------------------------------

      jrowstart = js + joff
      jrowend   = min(je + joff + 2,jmt)

!     for momentum... set coefficients for all latitudes

      if (first) then
!----------------------------------------------------------------------
!     anisotropic viscosity scheme of Large et al., 2001, JPO
!
!     coded by Christopher Somes
!     see Somes et al., 2010, GBC, auxiliary materials for additional
!     details
!----------------------------------------------------------------------
        istrt = 2
        iend = imt-1
        jstrt = 2
        jend = jmt-1

        V0 = 100.
        efD = 150000.
        aeddy = 1.e7
        N = 3.
        do jrow=jrowstart,jrowend
           do k=1,km
              do i=istrt,iend
                 beta = 0.0228e-11*abs(cos(pi/180*ulat(i,jrow)))
                 delx = dxudeg(i)*1.11e7
     &                  *abs(cos(pi/180*ulat(i,jrow)))
                 if((ulat(i,jrow).ge.-20).and.(ulat(i,jrow).le.20).and.
     &                (zw(k).le.55000)) then   ! apply only in tropics and upper ocean
                     if (umsk(i-1,jrow) == 0) then
                       WBC = 1.
                    else if (umsk(i-2,jrow) == 0)  then
                       WBC = 2.
                    else if (umsk(i-3,jrow) == 0) then
                       WBC = 3.
                    else if (umsk(i-4,jrow) == 0) then
                       WBC = 4.
                    else if (umsk(i-5,jrow) == 0) then
                       WBC = 5.
                    else if (umsk(i-6,jrow) == 0) then
                       WBC = 6.
                    else if (umsk(i-7,jrow) == 0) then
                       WBC = 7.
                    else if (umsk(i-8,jrow) == 0) then
                       WBC = 8.
                    else if (umsk(i-9,jrow) == 0) then
                       WBC = 9.
                    else if (umsk(i-10,jrow) == 0) then
                       WBC = 10.
                    else
                       WBC = 11.
                    endif
                    px = max(0., WBC - N)*delx/100000000.

                    bmunk = 0.2*beta*delx**3*exp(-px**2)
                    beddy = aeddy*(1+24.5*(1
     &                     - abs(cos(2*pi/180*ulat(i,jrow)))))

                    visc_cnu(i,k,jrow) = max(bmunk, beddy)
                 else
                    visc_cnu(i,k,jrow) = am
                 endif
              enddo
           enddo
        enddo
        do j=jstrt,jend
           jrow = j + joff
           do k=1,km
              do i=istrt-1,iend
                 gridlen = max(dxudeg(i)*1.11e7
     &                         *abs(cos(pi/180*ulat(i,j)))
     &                            , dyudeg(i)*1.1e7)
                 if ((ulat(i,j).ge.-20).and.(ulat(i,j).le.20).and.
     &                (zw(k).le.55000)) then  ! apply only in tropics and upper ocean
                    Vz = V0 !*exp(-1.*zw(k)/efD)
                    visc_ceu(i,k,j) = 0.5*Vz*gridlen
                 else
                    visc_ceu(i,k,j) = am
                 endif
              enddo
           enddo
        enddo
        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          do k=1,km
             do i=istrt,iend
                amc_north(i,k,jrow) = visc_cnu(i,k,jrow)*cst(jp1)
     &               *dytr(jp1)*csur(jrow)*dyur(jrow)
                amc_south(i,k,jrow) = visc_cnu(i,k,jrow)*cst(jrow)
     &               *dytr(jrow)*csur(jrow)*dyur(jrow)
             enddo
          enddo
        enddo
      endif

!     for tracers... set coefficients for all latitudes

      if (first) then
        diff_cnt  = ah
        diff_cet  = ah

        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          ahc_north(jrow) = diff_cnt*csu(jrow)*dyur(jrow)*cstr(jrow)
     &                              *dytr(jrow)
          ahc_south(jrow) = diff_cnt*csu(jm1)*dyur(jm1)*cstr(jrow)
     &                              *dytr(jrow)
        enddo
      endif

      return
      end
