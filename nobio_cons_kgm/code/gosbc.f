! source file: /raid24/aschmitt/UVic2.9/MOBI1.9/nobio/updates/gosbc.F
      subroutine gosbc (is, ie, js, je)

!=======================================================================
!     calculate the average fluxes for next ocean time step
!=======================================================================

      implicit none

      integer ie, is, je, js, i, j, nc
      integer iem1, isp1, jem1, jsp1, k

      real f1, f1a, f1l, fh, fs, fwcflx, fwaflx, time
      real area, tarea, tsflx, rsocn, tmp, fg, fgs, sulphfac

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "coord.h"
      include "grdvar.h"
      include "tmngr.h"
      include "switch.h"
      include "cembm.h"
      include "atm.h"
      include "mw.h"
      include "ice.h"
      include "mtlm.h"
      include "levind.h"

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      f1 = 1./atatm
      fh = 2.389e-8/atatm
      fs = -socn/atatm

!-----------------------------------------------------------------------
!     calculate average net fluxes. convert heat flux to cal/cm**2/s
!     from kW/m**2 and fresh water flux (cm/s) to an apparent salt
!     flux (g/cm**2/s) using global ocean average salinity, socn
!-----------------------------------------------------------------------

      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
            sbc(i,j,ihflx) = sbc(i,j,ihflx) + fh*flux(i,j,isat)
!           add virtual fluxes of salinity
            sbc(i,j,isflx) = sbc(i,j,isflx) + fs*flux(i,j,ishum)

          else
            sbc(i,j,ihflx) = 0.
            sbc(i,j,isflx) = 0.
          endif
          if (umsk(i,j) .ge. 0.5) then
            sbc(i,j,itaux) = f1*flux(i,j,nat+1)
            sbc(i,j,itauy) = f1*flux(i,j,nat+2)
          else
            sbc(i,j,itaux) = 0.
            sbc(i,j,itauy) = 0.
          endif
        enddo
      enddo

      call setbcx (sbc(1,1,ihflx), imt, jmt)
      call setbcx (sbc(1,1,isflx), imt, jmt)
      call setbcx (sbc(1,1,itaux), imt, jmt)
      call setbcx (sbc(1,1,itauy), imt, jmt)

!-----------------------------------------------------------------------
!     update boundary conditions from the land model
!     do this now instead of in gasbc so fields can be written out
!-----------------------------------------------------------------------

      f1l = 0.
      f1a = 0.
      if (atatm .ne. 0.) f1a = 1.0/atatm
      if (atlnd .ne. 0.) f1l = 1.0/atlnd
      do j=2,jmtm1
        do i=2,imtm1
          if (land_map(i,j) .ne. 0) then
            sbc(i,j,iro) = sbc(i,j,iro)*f1l
            sbc(i,j,isca) = sbc(i,j,isca)*f1l
            sbc(i,j,ievap) = sbc(i,j,ievap)*f1l
            sbc(i,j,ilwr) = sbc(i,j,ilwr)*f1l
            sbc(i,j,isens) = sbc(i,j,isens)*f1l
          else
            sbc(i,j,iro) = sbc(i,j,iro)*f1a
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,isca), imt, jmt)
      call setbcx (sbc(1,1,ievap), imt, jmt)
      call setbcx (sbc(1,1,ilwr), imt, jmt)
      call setbcx (sbc(1,1,isens), imt, jmt)
      call setbcx (sbc(1,1,iro), imt, jmt)

!-----------------------------------------------------------------------
!     zero diagnostic for river discharge and call river model
!-----------------------------------------------------------------------
      disch(:,:) = 0.
      call rivmodel

      return
      end

