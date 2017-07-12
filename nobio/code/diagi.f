! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/diagi.F
      subroutine diagi

!=======================================================================
!     initialize diagnostics quantities
!=======================================================================

      implicit none

      integer iobadt, iobads, jrow, k, ll, n, nreg, m, mask

      real zmau, zmat, zmsmf, zmsm, zmstf, zmst

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "stab.h"
      include "ctmb.h"
      include "coord.h"
      include "ctavg.h"
      include "diag.h"
      include "iounit.h"
      include "switch.h"
      include "tmngr.h"

      if (stabts) then
        numcfl = 0
        cflup  = c0
        cflvp  = c0
        cflwtp = c0
        cflwup = c0
        reynx  = c0
        reyny  = c0
        reynz  = c0
        peclx  = c0
        pecly  = c0
        peclz  = c0
        call getunit (iostab, 'iostab', 'formatted sequential rewind')
        rewind iostab
        endfile iostab
        rewind iostab
        call relunit (iostab)
        call getunit (iobadt, 'iobadt', 'formatted sequential rewind')
        rewind iobadt
        endfile iobadt
        rewind iobadt
        call relunit (iobadt)
        call getunit (iobads, 'iobads', 'formatted sequential rewind')
        rewind iobads
        endfile iobads
        rewind iobads
        call relunit (iobads)
      endif

      if (glents) then
        do jrow=1,jmt
          do k=1,km
            wtlev(k,jrow) = c0
            wulev(k,jrow) = c0
          enddo
          do k=0,km
            buoy(k,jrow) = c0
          enddo
          do ll=1,8
            do k=0,km
              engint(k,ll,jrow) = c0
            enddo
          enddo
          do ll=1,8
            engext(ll,jrow) = c0
          enddo
          tcerr(jrow) = c0
          ucerr(jrow) = c0
          itcerr(jrow) = 1
          jtcerr(jrow) = 1
          ktcerr(jrow) = 1
          iucerr(jrow) = 1
          jucerr(jrow) = 1
          kucerr(jrow) = 1
          wtbot(jrow)  = c0
          wubot(jrow)  = c0
          iwtbot(jrow) = 1
          jwtbot(jrow) = 1
          kwtbot(jrow) = 1
          iwubot(jrow) = 1
          jwubot(jrow) = 1
          kwubot(jrow) = 1
        enddo
      endif

      if (trmbts) then
        do n=1,nt
          ustf(n,1) = ' unknown  '
          ustf(n,2) = ' unknown units '
          if (n .eq. 1) then
            ustf(n,1) = ' stf(1)   = '
            ustf(n,2) = ' cal/cm**2/sec '
          endif
          if (n .eq. 2) then
            ustf(n,1) = ' stf(2)   = '
            ustf(n,2) = '     cm/sec '
          endif
        enddo

        do nreg=0,numreg
          if (nreg .gt. 0) then
            avgw(nreg) = c0
            do ll=1,17
              do k=0,km
                termbm(k,ll,1,nreg) = c0
                termbm(k,ll,2,nreg) = c0
              enddo
            enddo
          endif
          do n=1,nt
            do ll=1,15
              do k=0,km
                termbt(k,ll,n,nreg) = c0
              enddo
            enddo
          enddo
        enddo

        do nreg=0,nhreg
          smflx(1,nreg) = c0
          smflx(2,nreg) = c0
          do n=1,nt
            stflx(n,nreg) = c0
            asst(n,nreg)  = c0
          enddo
        enddo
      endif

      if (gyrets) then
        do jrow=1,jmt
          do m=1,ntmin2
            do ll=1,8
              ttn(ll,jrow,m) = c0
            enddo
          enddo
        enddo

        do n=0,nhreg
          do m=1,nt
            do jrow=1,jmt
              do ll=6,8
                ttn2(ll,jrow,m,n) = c0
              enddo
            enddo
          enddo
        enddo
      endif

      if (vmsfts) then
        do jrow=1,jmt
          do k=1,km
            vmsf(jrow,k) = c0
          enddo
        enddo
      endif

      if (tavgts) then
        do n=1,nt
          do mask=1,nhreg
            sumbf(mask,n) = c0
            do k=1,km
              sumbk(mask,k,n) = c0
              avgbk(mask,k,n) = c0
            enddo
          enddo
        enddo
      endif

      if (tsiperts .and. eots) then
        do jrow=1,jmt
          do k=0,km
            ektot(k,jrow) = c0
            do n=1,nt
              dtabs(k,n,jrow)  = c0
              tbar(k,n,jrow)   = c0
              travar(k,n,jrow) = c0
            enddo
          enddo
        enddo
      endif

      return
      end
