! source file: /usr/local/models/UVic_ESCM/2.9/source/common/filt.F
      subroutine filt (joff, js, je)

!=======================================================================
!     subroutine filt sets up input needed for fourier filtering
!     (when the "fourfil" option is defined) -or- symmetric finite
!     impulse response filtering (when the "firfil" option is defined)
!     of tracers at the specifiied high latitude row "jrow".
!=======================================================================

      implicit none

      integer n, j, js, je, jrow, joff, jj, isave, ieave, l, k
      integer is, ie, iredo, im, m, mm, idx, ism1, iea, i, ieb
      integer ii, jsf, jef

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "index.h"
      include "levind.h"
      include "mw.h"

      real tempik(imt,km)

      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!---------------------------------------------------------------------
!     fourier filter tracers at high latitudes
!---------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
      if ((jrow.gt.jft1.and.jrow.lt.jft2) .or. jrow.lt.jfrst) goto 101
      jj = jrow-jfrst+1
      if (jrow .ge. jft2) jj = jj-jskpt+1

!    if previous strips were of same length, do not recompute
!    fourier coeffs

      isave = 0
      ieave = 0
      do l=1,lsegf
        do k=1,km
          if (istf(jj,l,k) .ne. 0) then
            is    = istf(jj,l,k)
            ie    = ietf(jj,l,k)
            iredo = 0
            if (is.ne.isave .or. ie.ne.ieave) then
              iredo = -1
              isave = is
              ieave = ie
              im = ie-is+1
              if (im.ne.imtm2 .or. kmt(1,jrow).lt.k) then
                m = 1
                n = nint(im*cst(jrow)*cstr(jft0))
              else
                m = 3
                n = nint(im*cst(jrow)*cstr(jft0)*0.5)
              endif
            endif
            do mm=1,nt
              idx  = iredo+mm
              ism1 = is-1
              iea  = ie
              if (ie .ge. imt) iea = imtm1
              do i=is,iea
                tempik(i-ism1,k) = t(i,k,j,mm,taup1)
              enddo
              if (ie .ge. imt) then
                ieb = ie-imtm2
                ii  = imtm1-is
                do i=2,ieb
                  tempik(i+ii,k) = t(i,k,j,mm,taup1)
                enddo
              endif

              call filtr (tempik(1,k), im, m, n, idx)

              do i=is,iea
                t(i,k,j,mm,taup1) = tempik(i-ism1,k)
              enddo
              if (ie .ge. imt) then
                do i=2,ieb
                  t(i,k,j,mm,taup1) = tempik(i+ii,k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
101   continue
      enddo

      return
      end

