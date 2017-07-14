! source file: /raid24/aho/UVic2.9/default_comb2/nobio/updates/tracer_adv_flx.F
      subroutine adv_flux (joff, js, je, is, ie, n)

!=======================================================================
!     3rd order advective tracer flux

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!        n   = tracer
!=======================================================================

      implicit none

      integer istrt, iend, i, k, j, ip, kr, jq, lag, js, je, ip2, n
      integer joff, jrow, jp2, ie, is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "mw.h"

      parameter (istrt=2, iend=imt-1)

      real totvel, upos, uneg, vpos, vneg, wpos, wneg

      include "isopyc.h"

      lag = tau

!-----------------------------------------------------------------------
!     calculate 2*advective flux across eastern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km
          do i=istrt,iend-1
            ip2 = i + 2
            totvel = adv_vet(i,k,j)
     &              + adv_vetiso(i,k,j)
            upos = p5*(totvel + abs(totvel))
     &             *tmask(i-1,k,j)*tmask(i,k,j)*tmask(i+1,k,j)
            uneg = p5*(totvel - abs(totvel))
     &             *tmask(ip2,k,j)*tmask(i+1,k,j)*tmask(i,k,j)

            adv_fe(i,k,j) = totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau)
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)*t(i-1,k,j,n,lag))
     &                  - uneg*(curv_xn(i,1)*t(ip2,k,j,n,lag)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
          enddo
        enddo

        do k=1,km
          i=iend
          ip2 = 3
          totvel = adv_vet(i,k,j)
     &              + adv_vetiso(i,k,j)
          upos = p5*(totvel + abs(totvel))
     &             *tmask(i-1,k,j)*tmask(i,k,j)*tmask(i+1,k,j)
          uneg = p5*(totvel - abs(totvel))
     &             *tmask(ip2,k,j)*tmask(i+1,k,j)*tmask(i,k,j)

          adv_fe(i,k,j) = totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau)
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)*t(i-1,k,j,n,lag))
     &                  - uneg*(curv_xn(i,1)*t(ip2,k,j,n,lag)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
        enddo
        call setbcx (adv_fe(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across northern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      if (joff +js .eq. 2) then
        do j=1,1
          do k=1,km
            do i=2,imt-1
              adv_f4n(i,k,j,n) = c0
            enddo
          enddo
        enddo
      endif
      do j=js,je
        jrow = j + joff
        jp2 = min(j+2+joff,jmt) - joff
        do k=1,km
          do i=istrt,iend
            totvel = adv_vnt(i,k,j)
     &              + adv_vntiso(i,k,j)
            vpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k,j-1)*tmask(i,k,j)*tmask(i,k,j+1)
            vneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k,jp2)*tmask(i,k,j+1)*tmask(i,k,j)

            adv_f4n(i,k,j,n) = totvel*(
     &                          quick_y(jrow,1)*t(i,k,j  ,n,tau)
     &                        + quick_y(jrow,2)*t(i,k,j+1,n,tau))
     &                  - vpos*(curv_yp(jrow,1)*t(i,k,j+1,n,lag)
     &                         +curv_yp(jrow,2)*t(i,k,j  ,n,lag)
     &                         +curv_yp(jrow,3)*t(i,k,j-1,n,lag))
     &                  - vneg*(curv_yn(jrow,1)*t(i,k,jp2,n,lag)
     &                         +curv_yn(jrow,2)*t(i,k,j+1,n,lag)
     &                         +curv_yn(jrow,3)*t(i,k,j  ,n,lag))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across bottom face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=2,km-2
          do i=istrt,iend
            totvel = adv_vbt(i,k,j)
     &              + adv_vbtiso(i,k,j)
            wpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k+2,j)*tmask(i,k+1,j)*tmask(i,k,j)
            wneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k-1,j)*tmask(i,k,j)*tmask(i,k+1,j)

            adv_fb(i,k,j)  = totvel*(
     &                          quick_z(k,1)*t(i,k  ,j,n,tau)
     &                        + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                  - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                         +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                         +curv_zp(k,3)*t(i,k-1,j,n,lag))
     &                  - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                         +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                         +curv_zn(k,3)*t(i,k  ,j,n,lag))
          enddo
        enddo
        k=1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
     &            + adv_vbtiso(i,k,j)
          wpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k+2,j)*tmask(i,k+1,j)*tmask(i,k,j)

          adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau)
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                       +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                       +curv_zn(k,3)*t(i,k  ,j,n,lag))
        enddo
        k=km-1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
     &            + adv_vbtiso(i,k,j)
          wneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k-1,j)*tmask(i,k,j)*tmask(i,k+1,j)

          adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau)
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                       +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                       +curv_zp(k,3)*t(i,k-1,j,n,lag))
        enddo
      enddo

      return
      end
