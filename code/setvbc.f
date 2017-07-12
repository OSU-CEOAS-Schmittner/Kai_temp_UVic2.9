! source file: /raid24/aschmitt/UVic2.9/karin/mwc15_npzd_fe_n15_c13_alk_caco3/updates/setvbc.F
      subroutine setvbc (joff, js, je, is, ie)

!=======================================================================
!     set momentum and tracer vertical boundary conditions

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer js, je, istrt, is, iend, ie, n, j, i, jrow, joff, kz

      real uvmag

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "levind.h"
      include "scalar.h"
      include "mw.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

!----------------------------------------------------------------------
!       set no flux condition for all tracers at surface & bottom.
!----------------------------------------------------------------------

      do n=1,nt
        do j=js,je
          do i=istrt,iend
            stf(i,j,n) = c0
            btf(i,j,n) = c0
          enddo
        enddo
      enddo

!----------------------------------------------------------------------
!       apply surface tracer and momentum fluxes from the atmosphere
!       code is for 2 tracer and 2 momentum fluxes.
!----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,itemp) = sbc(i,jrow,ihflx)*tmask(i,1,j)
          stf(i,j,isalt) = sbc(i,jrow,isflx)*tmask(i,1,j)
          btf(i,j,itemp) = -bhf(i,jrow)*tmask(i,1,j)
          stf(i,j,idic) = sbc(i,jrow,idicflx)*tmask(i,1,j)
          stf(i,j,idic13) = sbc(i,jrow,idic13flx)*tmask(i,1,j)
          stf(i,j,ialk) = sbc(i,jrow,ialkflx)*tmask(i,1,j)
          stf(i,j,io2) = sbc(i,jrow,io2flx)*tmask(i,1,j)
          stf(i,j,ipo4) = sbc(i,jrow,ipo4flx)*tmask(i,1,j)
          stf(i,j,iphyt) = sbc(i,jrow,iphytflx)*tmask(i,1,j)
          stf(i,j,izoop) = sbc(i,jrow,izoopflx)*tmask(i,1,j)
          stf(i,j,idetr) = sbc(i,jrow,idetrflx)*tmask(i,1,j)
          stf(i,j,icocc) = sbc(i,jrow,icoccflx)*tmask(i,1,j)
          stf(i,j,icaco3) = sbc(i,jrow,icaco3flx)*tmask(i,1,j)
          stf(i,j,idop) = sbc(i,jrow,idopflx)*tmask(i,1,j)
          stf(i,j,ino3) = sbc(i,jrow,ino3flx)*tmask(i,1,j)
          stf(i,j,idon) = sbc(i,jrow,idonflx)*tmask(i,1,j)
          stf(i,j,idiaz) = sbc(i,jrow,idiazflx)*tmask(i,1,j)
          stf(i,j,idin15) = sbc(i,jrow,idin15flx)*tmask(i,1,j)
          stf(i,j,idon15) = sbc(i,jrow,idon15flx)*tmask(i,1,j)
          stf(i,j,iphytn15) = sbc(i,jrow,iphytn15flx)*tmask(i,1,j)
          stf(i,j,icoccn15) = sbc(i,jrow,icoccn15flx)*tmask(i,1,j)
          stf(i,j,izoopn15) = sbc(i,jrow,izoopn15flx)*tmask(i,1,j)
          stf(i,j,idetrn15) = sbc(i,jrow,idetrn15flx)*tmask(i,1,j)
          stf(i,j,idiazn15) = sbc(i,jrow,idiazn15flx)*tmask(i,1,j)
          stf(i,j,idfe) = sbc(i,jrow,idfeflx)*tmask(i,1,j)
          stf(i,j,idetrfe) = sbc(i,jrow,idetrfeflx)*tmask(i,1,j)
          stf(i,j,idoc13) = sbc(i,jrow,idoc13flx)*tmask(i,1,j)
          stf(i,j,iphytc13) = sbc(i,jrow,iphytc13flx)*tmask(i,1,j)
         stf(i,j,icoccc13) = sbc(i,jrow,icoccc13flx)*tmask(i,1,j)
         stf(i,j,icaco3c13) = sbc(i,jrow,icaco3c13flx)*tmask(i,1,j)
          stf(i,j,izoopc13) = sbc(i,jrow,izoopc13flx)*tmask(i,1,j)
          stf(i,j,idetrc13) = sbc(i,jrow,idetrc13flx)*tmask(i,1,j)
          stf(i,j,idiazc13) = sbc(i,jrow,idiazc13flx)*tmask(i,1,j)
          smf(i,j,1) = sbc(i,jrow,itaux)*umask(i,1,j)
          smf(i,j,2) = sbc(i,jrow,itauy)*umask(i,1,j)
        enddo
      enddo

!----------------------------------------------------------------------
!       set bottom drag
!----------------------------------------------------------------------

      do n=1,2
        if (cdbot .eq. c0) then
          do j=js,je
            do i=istrt,iend
              bmf(i,j,n) = c0
            enddo
          enddo
        else
          do j=js,je
            jrow = j + joff
            do i=istrt,iend
              kz = kmu(i,jrow)
              if (kz .ne. 0) then
                uvmag    = sqrt(u(i,kz,j,1,taum1)**2 +
     &                          u(i,kz,j,2,taum1)**2)
                bmf(i,j,n) = cdbot*u(i,kz,j,n,taum1)*uvmag
              else
                bmf(i,j,n) = c0
              endif
            enddo
          enddo
        endif
      enddo

!----------------------------------------------------------------------
!     apply zonal boundary conditions
!----------------------------------------------------------------------

      do n=1,nt
        call setbcx (stf(1,js,n), imt, je-js+1)
        call setbcx (btf(1,js,n), imt, je-js+1)
      enddo
      do n=1,2
        call setbcx (smf(1,js,n), imt, je-js+1)
        call setbcx (bmf(1,js,n), imt, je-js+1)
      enddo

      return
      end
