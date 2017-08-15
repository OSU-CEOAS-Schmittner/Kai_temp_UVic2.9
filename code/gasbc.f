! source file: /data/home/kai/dev/UVic2.9/updates/gasbc.F

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
c#  if defined 1
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
      atmpres = 1.0       !atm
      zero = 0.
!     fractionation factors from Zhang et al. 1995 Geochi.
!     Cosmochi. Acta
      ak = 0.99915 ! kinetic fractionation
      aaqg = 0.998764           ! aquatic - gas fractionation

      r13a = (dc13ccn*0.001 + 1.)*rc13std
      batmc13 = ak*aaqg*r13a

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
      sbc(:,:,idicflx) = 0.
      sbc(:,:,idic13flx) = 0.
      sbc(:,:,ialkflx) = 0.
      sbc(:,:,io2flx) = 0.
      sbc(:,:,ipo4flx) = 0.
      sbc(:,:,iphytflx) = 0.
      sbc(:,:,izoopflx) = 0.
      sbc(:,:,idetrflx) = 0.
      sbc(:,:,icoccflx) = 0.
      sbc(:,:,icaco3flx) = 0.
      sbc(:,:,idfeflx) = 0.
      sbc(:,:,idfeadep) = 0.
      sbc(:,:,idetrfeflx) = 0.
      sbc(:,:,idopflx) = 0.
      sbc(:,:,ino3flx) = 0.
      sbc(:,:,idonflx) = 0.
      sbc(:,:,idiazflx) = 0.
      sbc(:,:,idin15flx) = 0.
      sbc(:,:,idon15flx) = 0.
      sbc(:,:,iphytn15flx) = 0.
      sbc(:,:,icoccn15flx) = 0.
      sbc(:,:,izoopn15flx) = 0.
      sbc(:,:,idetrn15flx) = 0.
      sbc(:,:,idiazn15flx) = 0.
      sbc(:,:,iphytc13flx) = 0.
      sbc(:,:,icoccc13flx) = 0
      sbc(:,:,icaco3c13flx) = 0
      sbc(:,:,izoopc13flx) = 0.
      sbc(:,:,idetrc13flx) = 0.
      sbc(:,:,idoc13flx) = 0.
      sbc(:,:,idiazc13flx) = 0.

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
            sst = sbc(i,j,isst)
!           put reasonable limits on sst and sss for chemistry flux calculations
            sst = min(35.,max(sst,-2.))
            sss = min(45.,max(sss,0.))
            ao = 1. - aice(i,j,2)

!-----------------------------------------------------------------------
!           calculate ocean carbon fluxes
!-----------------------------------------------------------------------
            t_in = sst
            s_in = sss
            dic_in = sbc(i,j,issdic)
            ta_in = sbc(i,j,issalk)
            co2_in = co2ccn
            call co2calc_SWS (t_in, s_in, dic_in, ta_in, co2_in, atmpres
     &,                       zero, pH, co2star, dco2star, pCO2, dpco2
     &,                       CO3, Omega_c, Omega_a)
!           Schmidt number for CO2
            scco2 = 2073.1 - 125.62*sst + 3.6276*sst**2
     &            - 0.043219*sst**3
            piston_vel = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *((scco2/660.)**(-0.5))
!           dic in umol cm-3 or (mol m-3) => flux in umol cm-2 s-1
            sbc(i,j,idicflx) = piston_vel*dco2star
            adicg = 1.01051 - 1.05e-4*sst ! DIC-gas fractionation

            r13dic = sbc(i,j,issdic13)
     &               / (sbc(i,j,issdic)-sbc(i,j,issdic13))
            r13dic = min(r13dic, 2.*rc13std)
            r13dic = max(r13dic, 0.5*rc13std)
            bdic13 = ak*aaqg*r13dic/adicg

            sbc(i,j,idic13flx) = piston_vel
     &           *((batmc13/(1+batmc13))*(dco2star + co2star)
     &           - (bdic13/(1+bdic13))*co2star)
c            sbc(i,j,idic13flx) = sbc(i,j,idicflx)*rc13std/(1+rc13std)

!-----------------------------------------------------------------------
!           calculate ocean oxygen fluxes
!-----------------------------------------------------------------------
!           Schmidt number for O2
            sco2 = 1638.0 - 81.83*sst + 1.483*sst**2 - 0.008004*sst**3
!           piston velocity for O2
            piston_o2 = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *(sco2/660.0)**(-0.5)
!           oxygen saturation concentration [mol/m^3]
            f1 = alog((298.15 - sst)/(C2K + sst))
            f2 = f1*f1
            f3 = f2*f1
            f4 = f3*f1
            f5 = f4*f1
            o2sat = exp (2.00907 + 3.22014*f1 + 4.05010*f2
     &             + 4.94457*f3 - 2.56847E-1*f4 + 3.88767*f5
     &             + sss*(-6.24523e-3 - 7.37614e-3*f1 - 1.03410e-2*f2
     &             - 8.17083E-3*f3) - 4.88682E-7*sss*sss)
!           Convert from ml/l to mol/m^3
            o2sat = o2sat/22391.6*1000.0
            sbc(i,j,io2flx) = piston_o2*(o2sat - sbc(i,j,isso2))

          else
!-----------------------------------------------------------------------
!           calculate land carbon fluxes
!-----------------------------------------------------------------------
!           convert from kg m-2 s-1 => umol cm-2 s-1
            sbc(i,j,idicflx) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                       - sbc(i,j,iburn))*0.1/12.e-6
            sbc(i,j,idic13flx) = (sbc(i,j,inpp13) - sbc(i,j,isr13)
     &                       - sbc(i,j,iburn13))*0.1/12.e-6
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     set boundary conditions for carbon
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,idicflx), imt, jmt)
      call setbcx (sbc(1,1,idic13flx), imt, jmt)
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idicflx), dmsk, avgflxc)
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idic13flx), dmsk, avgflxc)
      carbemit = carbemit + co2emit*atmsa*segtim*daylen*1e-15

!-----------------------------------------------------------------------
!     set boundary conditions for oxygen
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,io2flx), imt, jmt)

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
