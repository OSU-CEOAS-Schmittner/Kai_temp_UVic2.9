! source file: /raid24/aschmitt/UVic2.9/MOBI1.9/nobio/updates/embmio.F
       subroutine embmout (is, ie, js, je)

!=======================================================================
!     output routine for energy-moisture balance model

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name
      character(32) :: nstamp
      character(3) :: a3

      integer i, ie, is, id_xt, id_yt, iou, j, je, js, n, L
      integer ndx, ntrec, it(10), ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical exists, inqvardef

      real avgper, ca, time, wt, wtp1
      real c100, c500, C2K, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "atm.h"
      include "solve.h"
      include "coord.h"
      include "grdvar.h"
      include "levind.h"
      include "ice.h"
      include "evp.h"
      include "mtlm.h"
      include "cembm.h"
      include "iounit.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "riv.h"
      include "cregin.h"

      real tmp_at(imt,jmt,nat)
      real sm(imt,jmt), st(imt,jmt), hs(imt,jmt), ro(imt,jmt)
      real rntatsl, rntatil
      real dmsk(imt,jmt)
      real tmp_dt(imt,jmt)

      real p_alb(is:ie,js:je), a_alb(is:ie,js:je), s_alb(is:ie,js:je)
      real sat(is:ie,js:je), ta_outswr(is:ie,js:je)
      real ta_netrad(is:ie,js:je)

      c100 = 100.
      c500 = 500.
      C2K = 273.15

      if (tsits .and. ntatia .ne. 0) then

!-----------------------------------------------------------------------
!     write atmospheric time series data
!-----------------------------------------------------------------------

        call ta_embm_tsi (is, ie, js, je, 2)

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (*,'(1x, a3, i7, 1x, a32, 3(a,1pe13.6))')
     &      'ts=',itt, nstamp, ' iterations =', tai_maxit
     &,     ' TAbar=', tai_sat, ' QAbar=', tai_shum
!begin AHO
          write (*,'(a,i5)') ' Index # [ntrec, embmio, ts=]: ', ntrec
!end AHO
        endif
        rntatil = 0.
        if (ntatil .gt. 0) rntatil = 1./float(ntatil)
        tai_hsno = tai_hsno + tai_LYING_SNOW*1000.*rntatil/rhosno
        call def_tsi
        call def_tsi_embm (fname)
        tai_catm = 0.

        avgper = tsiper*accel
        if (avgper .le. 1e-6) avgper = 0.
        tmp = 0.5
        time = time - tmp*avgper/365.
!       convert emissions from g cm-2 to kg s-1
        tmp = tai_co2emit*atmsa*1.e-3
!begin AHO
        write (*,'(a,i5)') 'Index # [ntrec, embmio], pre embm_tsi_out: 
     &    ', ntrec
!end AHO

        call embm_tsi_out (fname, avgper, time, nstamp, tai_sat
     &,                    tai_shum, tai_precip, tai_evap, tai_ohice
     &,                    tai_oaice, tai_hsno, tai_lhice, tai_laice
     &,                    tai_co2ccn, tmp, tai_dc14ccn
     &,                    tai_dc13ccn, tai_cfc11ccn
     &,                    tai_cfc12ccn, tai_maxit, tai_nsat, tai_ssat
     &,                    tai_nshum, tai_sshum, tai_nprecip
     &,                    tai_sprecip, tai_nevap, tai_sevap, tai_nohice
     &,                    tai_sohice, tai_noaice, tai_soaice, tai_nhsno
     &,                    tai_shsno, tai_nlhice, tai_slhice, tai_nlaice
     &,                    tai_slaice, tai_lsat, tai_osat, tai_lprecip
     &,                    tai_oprecip, tai_levap, tai_oevap
     &,                    tai_solins, tai_upsens, tai_uplwr, tai_outlwr
     &,                    tai_dnswr, tai_absswr, tai_netrad, tai_palb
     &,                    tai_aalb, tai_salb, tai_lsalb, tai_osalb
     &,                    tai_sst, tai_sss, tai_ssdic, tai_ssdic13
     &,                    tai_ssc14, tai_ssalk, tai_sso2, tai_sspo4
     &,                    tai_ssdop, tai_ssno3, tai_ssdon, tai_ssdin15
     &,                    tai_ssdon15, tai_ssdoc13, tai_sscfc11
     &,                    tai_sscfc12, tai_sulph, tai_volc, tai_agg
     &,                    tai_catm, tai_carbemit, ntrec
     &,                    tai_ssdfe)
!begin AHO
        write (*,'(a,i5)') 'Index # [ntrec, embmio], post embm_tsi_out:
     &    ', ntrec
!end AHO



        call ta_embm_tsi (is, ie, js, je, 0)

      endif

      if (timavgts .and. ntatsa .ne. 0) then

!-----------------------------------------------------------------------
!       write atmospheric time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_embm_tavg (is, ie, js, je, 2)

!       write time averaged data

        rntatsl = 0.
        if (ntatsl .gt. 0) rntatsl = 1./float(ntatsl)
        call unloadland (POINTS, TA_M, imt, jmt, land_map, sm)
        call unloadland (POINTS, TA_TSTAR_GB, imt, jmt, land_map, st)
        call unloadland (POINTS, TA_LYING_SNOW, imt, jmt, land_map, hs)
        call unloadland (POINTS, TA_SURF_ROFF, imt, jmt, land_map, ro)
        do j=js,je
          do i=is,ie
            sat(i,j) =  ta_at(i,j,isat) - elev(i,j)*rlapse
     &                - hicel(i,j,2)*rlapse
            if (ta_solins(i,j) .gt. 0) then
!             calculate incoming swr reaching the surface
              tmp = ta_solins(i,j) - ta_arswr(i,j) - ta_absin(i,j)
              s_alb(i,j) = 1. - ta_dnswr(i,j)/tmp
              a_alb(i,j) = ta_arswr(i,j)/ta_solins(i,j)
!             calculate outgoing swr leaving the atmosphere
              tmp = tmp - ta_dnswr(i,j) - ta_absout(i,j)
              ta_outswr(i,j) = ta_arswr(i,j) + tmp
              p_alb(i,j) = (ta_outswr(i,j))/ta_solins(i,j)
              ta_netrad(i,j) = ta_solins(i,j) - ta_outswr(i,j)
     &                       - ta_outlwr(i,j)
            else
              s_alb(i,j) = -1.
              a_alb(i,j) = -1.
              p_alb(i,j) = -1.
              ta_outswr(i,j) = 0.
              ta_netrad(i,j) = -ta_outlwr(i,j)
            endif
            if (land_map(i,j) .eq. 0) then
              sm(i,j) = ta_soilm(i,j)
              st(i,j) = ta_surf(i,j)
            else
!             convert from kg m-2 to cm (converted back later)
              sm(i,j) = sm(i,j)*rntatsl*0.1
!             convert from K to C (converted back later)
              st(i,j) = st(i,j)*rntatsl - C2K
!             convert from kg/m2/s to g/cm2/s (converted back later)
              ta_runoff(i,j) = ro(i,j)*rntatsl*0.1
      !       convert from kg/m2 to cm (converted back later)
              ta_hsno(i,j) = hs(i,j)*rntatsl*0.1/rhosno
            endif
            if (ta_hsno(i,j) .lt. 0) ta_hsno(i,j) = 0.
            if (ta_hice(i,j) .lt. 0) ta_hice(i,j) = 0.
            if (ta_aice(i,j) .lt. 0) ta_aice(i,j) = 0.
            if (ta_aice(i,j) .gt. 1) ta_aice(i,j) = 1.
          enddo
        enddo

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        call def_tavg
        call def_tavg_embm (fname)
        avgper = timavgper*accel
        if (avgper .le. 1e-6) avgper = 0.
        tmp = 0.5
        time = time - tmp*avgper/365.

        call embm_tavg_out (fname, is, ie, js, je, imt, jmt, nat, ncat
     &,   xt, yt, xu, yu, dxt, dyt, dxu, dyu, avgper, time, nstamp
     &,   mapat, ta_at, sat, ta_precip, ta_evap, ta_disch, ta_vflux
     &,   ta_outlwr, ta_uplwr, ta_upsens, ta_dnswr, ta_solins
     &,   ta_outswr, ta_netrad, p_alb, a_alb, s_alb, elev, ta_psno
     &,   ta_ws, ta_runoff
     &,   ta_wx, ta_wy
     &,   sm, st
     &,   ta_tice, ta_hice, ta_aice, ta_hsno
     &,   ta_uice, ta_vice, ta_xint, ta_yint
     &,   tlat, tlon, ulat, ulon, tgarea, ugarea
     &,   ta_dn, ta_de
     &,   ta_aicel, ta_hicel
     &,   tmsk, mskhr, nriv, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Atm time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp
!begin AHO
        write(*,'(a,i5)') ' Index # [ntrec, embmio, atm]: ]', ntrec
!end AHO


!       zero time average accumulators

        call ta_embm_tavg (is, ie, js, je, 0)

      endif

!-----------------------------------------------------------------------
!     do things that need to be done only once a year
!-----------------------------------------------------------------------

      if (eoyear) then

!       make new annual and running average
        call embmbc (atbar)
        if (totaltime .ge. yrlen*daylen - dtatm*0.5) then
          wt = 0.05
          do j=js,je
            do i=is,ie
              rtbar(i,j) = (1.-wt)*rtbar(i,j) + wt*atbar(i,j)/totaltime
              atbar(i,j) = 0.0
            enddo
          enddo
        endif
        totaltime = 0.0

        dmsk(:,:) = 1.
        tmp_dt(:,:) = rtbar(:,:) - tbar(:,:)
        call areaavg (tmp_dt, dmsk, dtbar)

      endif

!-----------------------------------------------------------------------
!     write atmospheric restart
!-----------------------------------------------------------------------

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_embm (0, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_embm (1, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
      endif

      return
      end

      subroutine ta_embm_tavg (is, ie, js, je, iflag)

!=======================================================================
!     atmospheric data time averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer is, ie, js, je, iflag, i, j, k, n, ndx

      real rntatsa

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "atm.h"
      include "riv.h"
      include "ice.h"
      include "csbc.h"

      real dmsk(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatsa = 0
        ta_at(:,:,:) = 0.
        ta_solins(:,:) = 0.
        ta_arswr(:,:) = 0.
        ta_dnswr(:,:) = 0.
        ta_absin(:,:) = 0.
        ta_absout(:,:) = 0.
        ta_uplwr(:,:) = 0.
        ta_upsens(:,:) = 0.
        ta_outlwr(:,:) = 0.
        ta_precip(:,:) = 0.
        ta_psno(:,:) = 0.
        ta_evap(:,:) = 0.
        ta_disch(:,:) = 0.
        ta_vflux(:,:) = 0.
        ta_ws(:,:) = 0.
        ta_runoff(:,:) = 0.
        ta_wx(:,:,:) = 0.
        ta_wy(:,:,:) = 0.
        ta_soilm(:,:) = 0.
        ta_surf(:,:) = 0.
        ta_tice(:,:) = 0.
        ta_hice(:,:) = 0.
        ta_aice(:,:) = 0.
        ta_hsno(:,:) = 0.
        ta_uice(:,:) = 0.
        ta_vice(:,:) = 0.
        ta_pice(:,:) = 0.
        ta_xint(:,:) = 0.
        ta_yint(:,:) = 0.
        ta_dn(:,:,:) = 0.
        ta_de(:,:,:) = 0.
        ta_aicel(:,:) = 0.
        ta_hicel(:,:) = 0.
        ta_psum(:) = 0.
      elseif (iflag .eq. 1) then

!       accumulate

        ntatsa = ntatsa + 1
        do n=1,nat
          ta_at(:,:,n) = ta_at(:,:,n) + at(:,:,2,n)
        enddo
        ta_solins(:,:) = ta_solins(:,:) + solins(:,:)
        ta_arswr(:,:) = ta_arswr(:,:) + solins(:,:)
     &                  *(1. - sbc(:,:,iaca))
        ta_dnswr(:,:) = ta_dnswr(:,:) + dnswr(:,:)
        ta_absin(:,:) = ta_absin(:,:) + solins(:,:)*sbc(:,:,iaca)
     &                  *scatter
        ta_absout(:,:) = ta_absout(:,:) + (solins(:,:)*sbc(:,:,iaca)
     &                   *pass - dnswr(:,:))*scatter
        ta_uplwr(:,:) = ta_uplwr(:,:) + uplwr(:,:)
        ta_upsens(:,:) = ta_upsens(:,:) + upsens(:,:)
        ta_outlwr(:,:) = ta_outlwr(:,:) + outlwr(:,:)
        ta_precip(:,:) = ta_precip(:,:) + precip(:,:)
        ta_psno(:,:) = ta_psno(:,:) + psno(:,:)
        ta_evap(:,:) = ta_evap(:,:) + evap(:,:)
        ta_disch(:,:) = ta_disch(:,:) + disch(:,:)
        ta_vflux(:,:) = ta_vflux(:,:) + vflux(:,:)
        ta_ws(:,:) = ta_ws(:,:) + sbc(:,:,iws)
        ta_runoff(:,:) = ta_runoff(:,:) + runoff(:,:)
        ta_wx(:,:,ishum) = ta_wx(:,:,ishum) + sbc(:,:,iwxq)
        ta_wy(:,:,ishum) = ta_wy(:,:,ishum) + sbc(:,:,iwyq)
        ta_wx(:,:,isat) = ta_wx(:,:,isat) + sbc(:,:,iwxt)
        ta_wy(:,:,isat) = ta_wy(:,:,isat) + sbc(:,:,iwyt)
        ta_soilm(:,:) = ta_soilm(:,:) + soilm(:,:,2)
        ta_surf(:,:) = ta_surf(:,:) + surf(:,:)
        ta_tice(:,:) = ta_tice(:,:) + tice(:,:)
        ta_hice(:,:) = ta_hice(:,:) + hice(:,:,2)
        ta_aice(:,:) = ta_aice(:,:) + aice(:,:,2)
        ta_hsno(:,:) = ta_hsno(:,:) + hsno(:,:,2)
        ta_uice(:,:) = ta_uice(:,:) + uice(:,:)
        ta_vice(:,:) = ta_vice(:,:) + vice(:,:)
        ta_pice(:,:) = ta_pice(:,:) + pice(:,:)
        ta_xint(:,:) = ta_xint(:,:) + xint(:,:)
        ta_yint(:,:) = ta_yint(:,:) + yint(:,:)
        ta_dn(:,:,:) = ta_dn(:,:,:) + dn(:,:,:)
        ta_de(:,:,:) = ta_de(:,:,:) + de(:,:,:)
        ta_aicel(:,:) = ta_aicel(:,:) + aicel(:,:,2)
        ta_hicel(:,:) = ta_hicel(:,:) + hicel(:,:,2)
        ta_psum(:) = ta_psum(:) + psum(:)

      elseif (iflag .eq. 2 .and. ntatsa .ne. 0) then

!       average
        rntatsa = 1./float(ntatsa)
        ta_at(:,:,:) = ta_at(:,:,:)*rntatsa
        ta_solins(:,:) = ta_solins(:,:)*rntatsa
        ta_arswr(:,:) = ta_arswr(:,:)*rntatsa
        ta_dnswr(:,:) = ta_dnswr(:,:)*rntatsa
        ta_absin(:,:) = ta_absin(:,:)*rntatsa
        ta_absout(:,:) = ta_absout(:,:)*rntatsa
        ta_uplwr(:,:) = ta_uplwr(:,:)*rntatsa
        ta_upsens(:,:) = ta_upsens(:,:)*rntatsa
        ta_outlwr(:,:) = ta_outlwr(:,:)*rntatsa
        ta_precip(:,:) = ta_precip(:,:)*rntatsa
        ta_psno(:,:) = ta_psno(:,:)*rntatsa
        ta_evap(:,:) = ta_evap(:,:)*rntatsa
        ta_disch(:,:) = ta_disch(:,:)*rntatsa
        ta_vflux(:,:) = ta_vflux(:,:)*rntatsa
        ta_ws(:,:) = ta_ws(:,:)*rntatsa
        ta_runoff(:,:) = ta_runoff(:,:)*rntatsa
        ta_wx(:,:,:) = ta_wx(:,:,:)*rntatsa
        ta_wy(:,:,:) = ta_wy(:,:,:)*rntatsa
        ta_soilm(:,:) = ta_soilm(:,:)*rntatsa
        ta_surf(:,:) = ta_surf(:,:)*rntatsa
        ta_hsno(:,:) = ta_hsno(:,:)*rntatsa
        ta_tice(:,:) = ta_tice(:,:)*rntatsa
        ta_hice(:,:) = ta_hice(:,:)*rntatsa
        ta_aice(:,:) = ta_aice(:,:)*rntatsa
        ta_uice(:,:) = ta_uice(:,:)*rntatsa
        ta_vice(:,:) = ta_vice(:,:)*rntatsa
        ta_pice(:,:) = ta_pice(:,:)*rntatsa
        ta_xint(:,:) = ta_xint(:,:)*rntatsa
        ta_yint(:,:) = ta_yint(:,:)*rntatsa
        ta_dn(:,:,:) = ta_dn(:,:,:)*rntatsa
        ta_de(:,:,:) = ta_de(:,:,:)*rntatsa
        ta_aicel(:,:) = ta_aicel(:,:)*rntatsa
        ta_hicel(:,:) = ta_hicel(:,:)*rntatsa
        ta_psum(:) = ta_psum(:)*rntatsa

      endif

      return
      end

      subroutine ta_embm_tsi (is, ie, js, je, iflag)

!=======================================================================
!     atmospheric data time integral averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer i, is, ie, j, js, je, iflag, maxit, n

      real rntatia, sum, tmp, tmp_solins, tmp_outlwr, tmp_dnswr

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "ice.h"
      include "solve.h"

      real sat(imt,jmt), dmsk(imt,jmt)
      real tmpij(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatia = 0
        tai_maxit = 0.
        tai_sat  = 0.
        tai_shum = 0.
        tai_precip = 0.
        tai_evap = 0.
        tai_co2ccn = 0.
        tai_co2emit = 0.
        tai_dc14ccn = 0.
        tai_dc13ccn = 0.
        tai_cfc11ccn = 0.
        tai_cfc12ccn = 0.
        tai_osat = 0.
        tai_oprecip = 0.
        tai_oevap = 0.
        tai_ohice = 0.
        tai_oaice = 0.
        tai_hsno = 0.
        tai_lsat = 0.
        tai_lprecip = 0.
        tai_levap = 0.
        tai_lhice = 0.
        tai_laice = 0.
        tai_nsat = 0.
        tai_nshum = 0.
        tai_nprecip = 0.
        tai_nevap = 0.
        tai_nohice = 0.
        tai_noaice = 0.
        tai_nhsno = 0.
        tai_nlhice = 0.
        tai_nlaice = 0.
        tai_ssat = 0.
        tai_sshum = 0.
        tai_sprecip = 0.
        tai_sevap = 0.
        tai_sohice = 0.
        tai_soaice = 0.
        tai_shsno = 0.
        tai_slhice = 0.
        tai_slaice = 0.
        tai_solins = 0.
        tai_dnswr = 0.
        tai_upsens = 0.
        tai_uplwr = 0.
        tai_outlwr = 0.
        tai_absswr = 0.
        tai_netrad = 0.
        tai_palb = 0.
        tai_aalb = 0.
        tai_salb = 0.
        tai_lsalb = 0.
        tai_osalb = 0.
        tai_sst = 0.
        tai_sss = 0.
        tai_ssdic = 0.
        tai_ssdic13 = 0.
        tai_ssc14 = 0.
        tai_ssalk = 0.
        tai_sso2 = 0.
        tai_sspo4 = 0.
        tai_ssdop = 0.
        tai_ssno3 = 0.
        tai_ssdon = 0.
        tai_ssdin15 = 0.
        tai_ssdon15 = 0.
        tai_ssdoc13 = 0.
        tai_sscfc11 = 0.
        tai_sscfc12 = 0.
        tai_sulph = 0.
        tai_volc = 0.
        tai_agg = 0.
        tai_carbemit = 0.
        tai_ssdfe = 0.

      elseif (iflag .eq. 1) then

!       accumulate
        ntatia = ntatia + 1
        maxit = 0
        do n=1,nat
          if (itout(n) .gt. maxit) maxit = itout(n)
        enddo
        tai_maxit = tai_maxit + maxit
        sat(:,:) = at(:,:,2,isat) - elev(:,:)*rlapse
     &           - hicel(:,:,2)*rlapse
!-----------------------------------------------------------------------
!       set data mask for global
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.

        call areaavg (sat, dmsk, tmp)
        tai_sat  = tai_sat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_shum = tai_shum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_precip = tai_precip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_evap = tai_evap + tmp
        tai_co2ccn = tai_co2ccn + co2ccn
        tai_co2emit = tai_co2emit + co2emit

!-----------------------------------------------------------------------
!       set data mask for global ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_osat  = tai_osat + tmp
        call areaavg (precip, dmsk, tmp)
        tai_oprecip  = tai_oprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_oevap  = tai_oevap + tmp
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_ohice = tai_ohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_oaice = tai_oaice + tmp
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_hsno = tai_hsno + tmp

!-----------------------------------------------------------------------
!       set data mask for global land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_lsat  = tai_lsat + tmp
        call areaavg (precip, dmsk, tmp)
        tai_lprecip  = tai_lprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_levap  = tai_levap + tmp
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_hsno = tai_hsno + tmp
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_lhice = tai_lhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_laice = tai_laice + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0.) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_nsat  = tai_nsat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_nshum  = tai_nshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_nprecip  = tai_nprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_nevap  = tai_nevap + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0. .or. tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_nohice = tai_nohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_noaice = tai_noaice + tmp
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_nhsno = tai_nhsno + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_nhsno = tai_nhsno + tmp
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_nlhice = tai_nlhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_nlaice = tai_nlaice + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0.) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_ssat  = tai_ssat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_sshum  = tai_sshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_sprecip  = tai_sprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_sevap  = tai_sevap + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0. .or. tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_sohice = tai_sohice + tmp
        call areatot (aice(1,1,2), dmsk, tmp)
        tai_soaice = tai_soaice + tmp
        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_shsno = tai_shsno + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areatot (hsno(1,1,2), dmsk, tmp)
        tai_shsno = tai_shsno + tmp
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areatot (hicel(1,1,2), dmsk, tmp)
        tai_slhice = tai_slhice + tmp
        call areatot (aicel(1,1,2), dmsk, tmp)
        tai_slaice = tai_slaice + tmp
        dmsk(:,:) = 1.
        call areaavg (solins, dmsk, tmp)
        tai_solins = tai_solins + tmp
!       save total incoming shortwave
        tmp_solins = tmp
        call areaavg (upsens, dmsk, tmp)
        tai_upsens = tai_upsens + tmp
        call areaavg (uplwr, dmsk, tmp)
        tai_uplwr = tai_uplwr + tmp
        call areaavg (outlwr, dmsk, tmp)
        tai_outlwr = tai_outlwr + tmp
!       save total outgoing longwave
        tmp_outlwr = tmp
        call areaavg (dnswr, dmsk, tmp)
        tai_dnswr = tai_dnswr + tmp
!       save total absorbed surface shortwave
        tmp_dnswr = tmp
!       shortwave not reflected by the atmosphere
        tmpij(:,:) = solins(:,:)*sbc(:,:,iaca)
        call areaavg (tmpij, dmsk, tmp)
!       atmospheric albedo is 1 - not refected/incoming
        if (tmp_solins .gt. 0) tai_aalb = tai_aalb+(1.-tmp/tmp_solins)
!       surface albedo is 1 - not refected/incoming
        if (tmp .gt. 0) tai_salb = tai_salb + (1.-tmp_dnswr/(tmp*pass))
!       total absobed shortwave for planetary albedo
!       this is confusing, it really is: dnswr +  solins*sbc(iaca)*scatter
!       + (solins*sbc(iaca)*pass - dnswr)*scatter
        tmpij(:,:) = tmpij(:,:)*scatter*(1. + pass)+dnswr(:,:)*pass
        call areaavg (tmpij, dmsk, tmp)
        tai_absswr = tai_absswr + tmp
!       planetary albedo is 1 - not refected/incoming
        if (tmp_solins .gt. 0) tai_palb = tai_palb+(1.-tmp/tmp_solins)
        tai_netrad = tai_netrad + tmp - tmp_outlwr
        dmsk(:,:) = 1.
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.
        call areaavg (dnswr, dmsk, tmp_dnswr)
        tmpij(:,:) = solins(:,:)*sbc(:,:,iaca)
        call areaavg (tmpij, dmsk, tmp)
        if (tmp .gt. 0) tai_lsalb = tai_lsalb+(1.-tmp_dnswr/(tmp*pass))
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areaavg (dnswr, dmsk, tmp_dnswr)
        call areaavg (tmpij, dmsk, tmp)
        if (tmp .gt. 0) tai_osalb = tai_osalb+(1.-tmp_dnswr/(tmp*pass))

        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areaavg (sbc(1,1,isst), dmsk, tmp)
        tai_sst = tai_sst + tmp
        call areaavg (sbc(1,1,isss), dmsk, tmp)
        tai_sss = tai_sss + tmp
        tai_carbemit = tai_carbemit + carbemit

      elseif (iflag .eq. 2 .and. ntatia .ne. 0) then

!       average
        rntatia = 0.
        if (ntatia .gt. 0.) rntatia = 1./float(ntatia)
        tai_maxit = tai_maxit*rntatia
        tai_sat  = tai_sat*rntatia
        tai_shum = tai_shum*rntatia
        tai_precip = tai_precip*rntatia
        tai_evap = tai_evap*rntatia
        tai_co2ccn = tai_co2ccn*rntatia
        tai_co2emit = tai_co2emit*rntatia
        tai_dc13ccn = tai_dc13ccn*rntatia
        tai_dc14ccn = tai_dc14ccn*rntatia
        tai_cfc11ccn = tai_cfc11ccn*rntatia
        tai_cfc12ccn = tai_cfc12ccn*rntatia
        tai_osat = tai_osat*rntatia
        tai_oprecip = tai_oprecip*rntatia
        tai_oevap = tai_oevap*rntatia
        tai_ohice = tai_ohice*rntatia
        tai_oaice = tai_oaice*rntatia
        tai_hsno = tai_hsno*rntatia
        tai_lsat = tai_lsat*rntatia
        tai_lprecip = tai_lprecip*rntatia
        tai_levap = tai_levap*rntatia
        tai_lhice = tai_lhice*rntatia
        tai_laice = tai_laice*rntatia
        tai_nsat = tai_nsat*rntatia
        tai_nshum = tai_nshum*rntatia
        tai_nprecip = tai_nprecip*rntatia
        tai_nevap = tai_nevap*rntatia
        tai_nohice = tai_nohice*rntatia
        tai_noaice = tai_noaice*rntatia
        tai_nhsno = tai_nhsno*rntatia
        tai_nlhice = tai_nlhice*rntatia
        tai_nlaice = tai_nlaice*rntatia
        tai_ssat = tai_ssat*rntatia
        tai_sshum = tai_sshum*rntatia
        tai_sprecip = tai_sprecip*rntatia
        tai_sevap = tai_sevap*rntatia
        tai_sohice = tai_sohice*rntatia
        tai_soaice = tai_soaice*rntatia
        tai_shsno = tai_shsno*rntatia
        tai_slhice = tai_slhice*rntatia
        tai_slaice = tai_slaice*rntatia
        tai_solins = tai_solins*rntatia
        tai_upsens = tai_upsens*rntatia
        tai_uplwr = tai_uplwr*rntatia
        tai_outlwr = tai_outlwr*rntatia
        tai_dnswr = tai_dnswr*rntatia
        tai_absswr = tai_absswr*rntatia
        tai_netrad = tai_netrad*rntatia
        tai_palb = tai_palb*rntatia
        tai_aalb = tai_aalb*rntatia
        tai_salb = tai_salb*rntatia
        tai_lsalb = tai_lsalb*rntatia
        tai_osalb = tai_osalb*rntatia
        tai_sst = tai_sst*rntatia
        tai_sss = tai_sss*rntatia
        tai_ssdic = tai_ssdic*rntatia
        tai_ssdic13 = tai_ssdic13*rntatia
        tai_ssc14 = tai_ssc14*rntatia
        tai_ssalk = tai_ssalk*rntatia
        tai_sso2 = tai_sso2*rntatia
        tai_sspo4 = tai_sspo4*rntatia
        tai_ssdop = tai_ssdop*rntatia
        tai_ssno3 = tai_ssno3*rntatia
        tai_ssdon = tai_ssdon*rntatia
        tai_ssdin15 = tai_ssdin15*rntatia
        tai_ssdon15 = tai_ssdon15*rntatia
        tai_ssdoc13 = tai_ssdoc13*rntatia
        tai_sscfc11 = tai_sscfc11*rntatia
        tai_sscfc12 = tai_sscfc12*rntatia
        tai_sulph = tai_sulph*rntatia
        tai_volc = tai_volc*rntatia
        tai_agg = tai_agg*rntatia
        tai_carbemit = tai_carbemit*rntatia
        tai_ssdfe = tai_ssdfe*rntatia

      endif

      return
      end
