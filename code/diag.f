! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/diag.F
      subroutine diag (joff, js, je, is, ie)

!=======================================================================
!     calculate diagnostics

!     input:

!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!      is     = starting longitude index for calculations
!      ie     = ending longitude index for calculations
!=======================================================================

      implicit none

      character(120) :: fname
      character(32) :: nstamp

      integer ntrec, nyear, nmonth, nday, nhour, nmin, nsec
      integer i, k, j, ip, kr, jq, js, je, istrt, is, iend, ie, joff
      integer jrow, n, jlat, jj, indp, ks, ke, m, io, i1, i2, iocm

      real zmau, zmat, zma1, zmsmf, zmsm, zma2, zmstf, zmst, reltim
      real fx, scl, period, ce, cn, cb, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "diaga.h"
      include "docnam.h"
      include "grdvar.h"
      include "iounit.h"
      include "isopyc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "vmixc.h"
      include "levind.h"
      include "emode.h"
      include "cembm.h"

      real tmp_t(imt,km,nt), tmp_stf(imt,nt)
      real vbarx(km)
      real aibuf(imt,km)

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit longitudes
!-----------------------------------------------------------------------

      istrt  = max(2,is)
      iend   = min(imt-1,ie)

      if (tsiperts .and. .not. euler2 .and. joff .eq. 0)
     &  nv_otsf = nv_otsf + 1
      if (tsiperts .and. .not. euler2 .and. joff .eq. 0)
     &  nt_slh = nt_slh + 1
      do j=js,je
        jrow = joff + j
!-----------------------------------------------------------------------
!       diagnostic: accumulate "tau" data for time step integrals
!-----------------------------------------------------------------------

        if (tsiperts .and. .not. euler2) then

          if (jrow .ge. jsot .and. jrow .le. jeot) then
            if (mrot .gt. 0 .and. mrot .le. nhreg) then
              do i=2,imtm1
                if (mskhr(i,jrow) .eq. mrot) then
                  do k=1,kmu(i,jrow)
                    v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                               dxu(i)
                  enddo
                endif
              enddo
            else
              do i=isot1,ieot1
                do k=1,kmu(i,jrow)
                  v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                             dxu(i)
                enddo
              enddo
              do i=isot2,ieot2
                do k=1,kmu(i,jrow)
                  v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                             dxu(i)
                enddo
              enddo
            endif
          endif
          do i=1,imt
            do k=1,kmt(i,jrow)
              t_slh(i,jrow,k,1) = t_slh(i,jrow,k,1) + t(i,k,j,1,tau)
              t_slh(i,jrow,k,2) = t_slh(i,jrow,k,2) + t(i,k,j,2,tau)
            enddo
          enddo

        endif

!-----------------------------------------------------------------------
!       diagnostic: accumulate "tau" data for time means
!-----------------------------------------------------------------------

        if (timavgperts .and. .not. euler2) then
          if (istrt .ne. 2 .and. iend .ne. imt-1) then
            write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "avgvar"'
            stop '=>diag'
          else
            call avgvar (j, jrow, adv_vbt(1,1,j), u(1,1,1,1,tau)
     &,                  t(1,1,1,1,tau), stf, smf, mapt)
          endif
        endif

!-----------------------------------------------------------------------
!       diagnostic: compute stability diagnostics
!-----------------------------------------------------------------------

        if (stabts .and. eots) then
          if (istrt .ne. 2 .and. iend .ne. imt-1) then
            write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "stab"'
            stop '=>diag'
          else
            call stab (j, jrow)
          endif
        endif

!-----------------------------------------------------------------------
!       construct meridional overturning of mass
!-----------------------------------------------------------------------

        if (jrow .lt. jmtm1 .and. vmsfts .and. eots) then
          do k=1,km
            vbarx(k) = c0
          enddo

          do k=1,km
            do i=istrt,iend
              vbarx(k) = vbarx(k) + u(i,k,j,2,tau)*csu(jrow)*dxu(i)
            enddo
            if (k .eq. 1) then
              vmsf(jrow,k) = vbarx(k)*dzt(k)
            else
              vmsf(jrow,k) = vmsf(jrow,k-1) + vbarx(k)*dzt(k)
            endif
          enddo
        endif
      enddo

      return
      end
