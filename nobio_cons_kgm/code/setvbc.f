! source file: /raid24/aho/UVic2.9/default_comb2/nobio/updates/setvbc.F
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
