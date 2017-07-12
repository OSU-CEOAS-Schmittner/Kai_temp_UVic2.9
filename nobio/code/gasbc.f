! source file: /raid24/aschmitt/UVic2.9/MOBI1.9/nobio/updates/gasbc.F

      subroutine gasbc (is, ie, js, je)

!=======================================================================
!     calculate boundary conditions for the atmospheric model
!=======================================================================

      implicit none

      integer ie, is, je, js, i, iem1, isp1, j, jem1, jsp1, k, n

      real sss, sst, xconv, t_in, s_in, dic_in, ta_in, co2_in
      real atmpres, pH, co2star,  dco2star, pCO2
      real dpco2, CO3, Omega_c, Omega_a, scco2, piston_vel, avgflxc
      real calday, f, sco2, o2sat, o2sato, o2surf, piston_o2, cfc11ccn
      real cfc12ccn, wt, sccfc, piston_cfc, sol_cfc, cfcsat, ao, tarea
      real tdc14ccn, h_r, d, f1, f2, f3, f4, f5, area, C2K, tmp, zero
      real depth
      real ak, aaqg, adicg, r13a, r13dic, batmc13, bdic13

      include "size.h"
      include "npzd.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "mw.h"
      include "ice.h"
      include "switch.h"
      include "tmngr.h"
      include "cembm.h"
      include "atm.h"
      include "insolation.h"
      include "calendar.h"
      include "grdvar.h"
      include "levind.h"
      include "solve.h"
      include "mtlm.h"
c#  if defined O_mtlm_carbon_13
c      include "mtlmc13.h"
c#  endif
c#  if defined O_mtlm_carbon_14
c      include "mtlmc14.h"
c#  endif

      real cosz(is:ie,js:je)
      real dmsk(is:ie,js:je)

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!     xconv is constant to convert piston_vel from cm/hr -> cm/s
!     here it is 100.*a*xconv (100 => m to cm, a=0.337, xconv=1/3.6e+05)
      xconv = 33.7/3.6e+05
      xconv = xconv*0.75
      C2K = 273.15

!-----------------------------------------------------------------------
!     zero totals for new accumulation
!-----------------------------------------------------------------------
      atatm = 0.
      flux(:,:,:) = 0.
      sbc(:,:,ihflx) = 0.
      sbc(:,:,isflx) = 0.
      sbc(:,:,iro) = 0.
      sbc(:,:,iat) = 0.
      sbc(:,:,irh) = 0.
      sbc(:,:,ipr) = 0.
      sbc(:,:,ips) = 0.
      sbc(:,:,iaws) = 0.
      sbc(:,:,iswr) = 0.

!-----------------------------------------------------------------------
!     set solar constant
!-----------------------------------------------------------------------
      call solardata

!-----------------------------------------------------------------------
!     update insolation for the current day
!-----------------------------------------------------------------------
!     subroutine decl is expecting a 365.25 day year
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      i = (ie-is+1)*(je-js+1)
!     get average zenith angle
      call zenith (i, c0, daylen, daylen, tlat, tlon, sindec, cosz)
      solins(is:ie,js:je) = solarconst*eccf*cosz(is:ie,js:je)

!-----------------------------------------------------------------------
!     set co2 concentration or emissions by tracking average co2
!-----------------------------------------------------------------------
      call co2ccndata

!-----------------------------------------------------------------------
!     update any atmospheric data
!-----------------------------------------------------------------------
      call data (is, ie, js, je)

!-----------------------------------------------------------------------
!     calculate freezing point of sea water using UNESCO (1983)
!-----------------------------------------------------------------------

      do j=jsp1,jem1
        do i=isp1,iem1

          if (tmsk(i,j) .ge. 0.5) then

            sss = 1000.0*sbc(i,j,isss) + 35.0
            frzpt(i,j) = -.0575*sss + 1.71e-3*sss**1.5 - 2.155e-4*sss**2
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate CO2 forcing
!-----------------------------------------------------------------------
      call co2forc

!-----------------------------------------------------------------------
!     set flags to calculate new coefficients
!-----------------------------------------------------------------------
      newcoef(:,:) = .true.

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_embm_tavg (is, ie, js, je, 0)

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_embm_tsi (is, ie, js, je, 0)

      return
      end
