! source file: /data/home/kai/dev/UVic2.9/updates/setembm.F
      subroutine setembm (is, ie, js, je)

!=======================================================================
!     initialize the energy-moisture balance model
!=======================================================================

      implicit none

      character(120) :: fname, vname, new_file_name, text
      character(3) :: a3

      integer i, ie, ii, iou, is, j, je, jj, js, jz, k, m, n, nsolve
      integer nu, nsum, ib(10), ic(10)

      logical exists, inqvardef

      real dlam, dphi, dtatms, dte, dyz, eccice, grarea, saltmax
      real si, ssh, t1, tair, yz_max, yz_min, wz, calday, tmp
      real zrel, c100, c1e4, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "solve.h"
      include "switch.h"
      include "coord.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "insolation.h"
      include "ice.h"
      include "evp.h"
      include "riv.h"
      include "tmngr.h"
      include "levind.h"
      include "csbc.h"
      include "scalar.h"
      include "veg.h"
      real rveg(imt,jmt)
      real dmsk(imt,jmt), tmpij(imtm2,jmtm2)
      real tmp_dt(imt,jmt)

      c100      = 100.
      c1e4      = 1.e4
      C2K       = 273.15

      cdatm     = 1.e-3
      cpatm     = 1.004e7
      sht       = 8.4e5
      shq       = 1.8e5
      shc       = 8.049e5
      rhoatm    = 1.250e-3
      esatm     = 4.6e-05
      pcfactor  = 0.
      cssh      = 3.8011e-3
      cfc11ccnn = 0.
      cfc11ccns = 0.
      cfc12ccnn = 0.
      cfc12ccns = 0.
      dc14ccnn  = 0.
      dc14ccne  = 0.
      dc14ccns  = 0.

      rhoocn    = 1.035
      esocn     = 5.4e-5
      vlocn     = 2.501e10

      cdice     = 5.5e-3
      rhoice    = 0.913
      rhosno    = 0.330
      esice     = 5.347e-5
      slice     = 2.835e10
      flice     = 3.34e9
      condice   = 2.1656e5

      soilmax   = 15.
      eslnd     = 5.347e-5

      nivc      = 1
      dtatms    = 1800.
      ns        = 30

      dalt_v    = 3.3e-3
      dalt_o    = 1.4e-3
      dalt_i    = 1.4e-3

!     ensure pass is between zero and one.
      pass =  min(max((1. - scatter), 0.), 1.)

!     gtoppm is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      gtoppm = 1./(4.138e-7*rhoatm*shc)

!     calculate atmospheric surface area
      atmsa = 0.
      do j=2,jmtm1
        do i=2,imtm1
          atmsa = atmsa + dxt(i)*dyt(j)*cst(j)
        enddo
      enddo

      if (mod(timavgint, segtim) .gt. 1.e-6 .and. timavgint .gt. 0.)
     &  then
        t1 = nint(timavgint/segtim)*segtim
        if (t1 .lt. segtim) t1 = segtim
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: "timavgint" does not contain an integral number'
     &,   '              of coupling time steps "segtim".              '
        write (stdout,*) '              (changed "timavgint" from '
     &, timavgint,' days to ', t1,' days to insure this condition)'
        timavgint = t1
      endif
      if (timavgint .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging interval "timavgint" = 0. implies no '
     &,   '             averaging when "time_averages" is enabled      '
      endif
      if (timavgint .gt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the interval "timavgint" exceeds the averaging '
     &,   '             period "timavgper" for option "time_averages"  '
      endif
      if (timavgint .lt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging period "timavgper" exceeds interval  '
     &,   '             "timavgint". Setting timavgper = timavgint     '
        timavgper = timavgint
      endif
      if (timavgper .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the averaging period "timavgper" is zero. The  '
     &,   '             average will be over only one time step!       '
      endif
      write (stdout,'(/,1x,a,f10.2,a,/,1x,a,f10.2,a)')
     &  '==> Time averages will be written every ', timavgint, ' days, '
     &, '    with an averaging period of         ', timavgper, ' days. '

!-----------------------------------------------------------------------
!     read CO2 concentration from data
!-----------------------------------------------------------------------
      call co2ccndata

!-----------------------------------------------------------------------
!     calculate the relative CO2 forcing term
!-----------------------------------------------------------------------

      call co2forc

      write (stdout,*)
      write (stdout,*) 'CO2 ratio (reference = 280 ppmv) =',co2ccn/280.
      write (stdout,*) 'Yields radiative forcing (W/m2) = ',anthro*1.e-3

!-----------------------------------------------------------------------
!     calculate the expansion coefficients for Berger's solution for
!     the year of the initial conditions
!-----------------------------------------------------------------------

      write (stdout,*)
      write (stdout,*) 'Initial Orbital Parameters:'
      call orbit (orbit_yr, eccen, obliq, mvelp, lambm0)
      write (stdout,*) '  Orbital Year:', orbit_yr
      write (stdout,*) '  Eccentricity:', eccen
      write (stdout,*) '  Obliquity:   ', obliq
      write (stdout,*) '  Longitude of Perihelion:', mvelp+180.

!-----------------------------------------------------------------------
!     calculate annual average insolation and Coriolis factor
!-----------------------------------------------------------------------

      radian = 360./(2.*pi)
      do j=1,jmt
        do i=1,imt
!         calculate coriolis parameter
          fcor(i,j) = 2.*omega*sin(ulat(i,j)/radian)
        enddo
      enddo

!-----------------------------------------------------------------------
!     read diffusion
!-----------------------------------------------------------------------

      dn(:,:,:) = 5.e9
      de(:,:,:) = 5.e9
      fname = new_file_name ("A_diff.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        call openfile (fname, iou)
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        do n=1,nat
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
!         northward component
          vname = 'dn_'//trim(a3)
          if (trim(mapat(n)) .eq. 'sat') then
            vname = 'A_difftY'
          elseif (trim(mapat(n)) .eq. 'shum') then
            vname = 'A_diffqY'
          elseif (trim(mapat(n)) .eq. 'co2') then
            vname = 'A_diffcY'
          endif
          exists = inqvardef(trim(vname), iou)
          if (exists) then
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, c1e4, c0)
            dn(2:imtm1,2:jmtm1,n) = tmpij(1:imtm2,1:jmtm2)
            call embmbc (dn(:,:,n))
          endif
!         eastward component
          vname = 'de_'//trim(a3)
          if (trim(mapat(n)) .eq. 'sat') then
            vname = 'A_difftX'
          elseif (trim(mapat(n)) .eq. 'shum') then
            vname = 'A_diffqX'
          elseif (trim(mapat(n)) .eq. 'co2') then
            vname = 'A_diffcX'
          endif
          exists = inqvardef(trim(vname), iou)
          if (exists) then
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, c1e4, c0)
            de(2:imtm1,2:jmtm1,n) = tmpij(1:imtm2,1:jmtm2)
            call embmbc (de(:,:,n))
          endif
        enddo
      endif

!-----------------------------------------------------------------------
!     set solver parameters
!-----------------------------------------------------------------------

      nsolve = 0
      itin(:)  = 500        ! max solver iterations
      epsin(:) = 5.e-7
      epsin(ishum) = 1.e-5
      epsin(isat) = 1.e-3
      nsolve = nsolve + 1
      levelin = 20              ! max coarse grid level
      if (nsolve .ne. 1) then
        write(*,*) '==> Error: more or less than one solver defined.'
        write(*,*) '           Use only one of embm_adi, embm_mgrid,'
     &,   ' embm_slap, embm_essl, embm_sparskit or embm_explicit'
        stop '=>setembm'
      endif

!-----------------------------------------------------------------------
!     check latent heats will sum to zero
!-----------------------------------------------------------------------

      if (slice .ne. vlocn + flice) write (stdout,'(/,a)')
     &   '==> Warning: changing latent heat of fusion to conserve heat'
        flice = slice - vlocn

!-----------------------------------------------------------------------
!     calculate grid terms for the atmospheric solver
!-----------------------------------------------------------------------

      do j=2,jmtm1
        dsgrd(j) = csu(j-1)/(dyu(j-1)*cst(j)*dyt(j))
        dngrd(j) = csu(j)/(dyu(j)*cst(j)*dyt(j))
        asgrd(j) = csu(j-1)/(2.*cst(j)*dyt(j))
        angrd(j) = csu(j)/(2.*cst(j)*dyt(j))
      enddo
      do i=2,imtm1
        dwgrd(i) = 1./(dxu(i-1)*dxt(i))
        degrd(i) = 1./(dxu(i)*dxt(i))
        azgrd(i) = 1./(2.*dxt(i))
      enddo

!-----------------------------------------------------------------------
!     set initial conditions or read a restart
!-----------------------------------------------------------------------

      newcoef(:,:) = .true.

      nats = namix
      dayoyr = 1.
      itt = 0
      irstdy = 0
      msrsdy = 0
      totaltime = 0.
      atbar(:,:) = 0.
      rtbar(:,:) = 0.
      at(:,:,:,:) = 0.
      tair = 13.
      at(:,:,:,isat) = tair
      ssh = cssh*exp(17.67*tair/(tair + 243.5))
      rh(:,:) = rhmax
      at(:,:,:,ishum) = rhmax*ssh
      carbemit = 0.
      precip(:,:) = 0.
      aicel(:,:,:) = 0.
      hicel(:,:,:) = 0.
      soilm(:,:,:) = 0.
      surf(:,:) = 0.
      hice(:,:,:) = 0.
      aice(:,:,:) = 0.
      tice(:,:) = 0.
      hsno(:,:,:) = 0.
      uice(:,:) = 0.
      vice(:,:) = 0.
      sbc(:,:,isu) = 0.
      sbc(:,:,isv) = 0.
      sbc(:,:,igu) = 0.
      sbc(:,:,igv) = 0.
      sig11n(:,:) = 0.
      sig11e(:,:) = 0.
      sig11s(:,:) = 0.
      sig11w(:,:) = 0.
      sig22n(:,:) = 0.
      sig22e(:,:) = 0.
      sig22s(:,:) = 0.
      sig22w(:,:) = 0.
      sig12n(:,:) = 0.
      sig12e(:,:) = 0.
      sig12s(:,:) = 0.
      sig12w(:,:) = 0.
      bv(:) = 0.
      xv(:) = 0.

!-----------------------------------------------------------------------
!     set land ice data and tracer grid ocean mask
!-----------------------------------------------------------------------

      call icedata
      dsealev = sealev

      if (.not. init) then
        fname = new_file_name ("restart_embm.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) call embm_rest_in (fname, is, ie, js, je)
      endif

!-----------------------------------------------------------------------
!     read average air temperature
!-----------------------------------------------------------------------

      tbar(:,:) = 0.
      fname = new_file_name ("A_slatref.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'A_slat in setembm.f'
      endif
      ib(:) = 1
      ic(:) = 1
      ic(1) = imtm2
      ic(2) = jmtm2
      call openfile (fname, iou)
      exists = inqvardef('A_slat', iou)
      if (.not. exists) then
        print*, "Error => A_slat does not exist."
        stop 'A_slat in setembm.f'
      else
        call getvara ('A_slat', iou, imtm2*jmtm2, ib, ic, tmpij, c1, c0)
        tbar(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        text = "C"
        call getatttext (iou, 'A_slat', 'units', text)
!       convert to model units (C)
        if (trim(text) .eq. "K")
     &    where (tbar(:,:) .lt. 1.e30) tbar(:,:) = tbar(:,:) - C2K
      endif
      call embmbc (tbar)

!-----------------------------------------------------------------------
!     read land elevations
!-----------------------------------------------------------------------

      elev(:,:) = 0.
      fname = new_file_name ("L_elev.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        call getvara ('L_elev', iou, imtm2*jmtm2, ib, ic, tmpij
     &,   c100, c0)
        elev(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (elev)
      endif
!     check for negative elevations
      where (elev(:,:) .lt. 0.) elev(:,:) = 0.

!-----------------------------------------------------------------------
!     initialize running annual averages
!-----------------------------------------------------------------------

      fname = new_file_name ("restart_embm.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('rtbar', iou)
      endif
      if (.not. exists .or. init) then
        totaltime = 0.
        atbar(:,:) = 0.
        rtbar(:,:) = tbar(:,:)
      endif

      dmsk(:,:) = 1.
      tmp_dt(:,:) = rtbar(:,:) - tbar(:,:)
      call areaavg (tmp_dt, dmsk, dtbar)

!-----------------------------------------------------------------------
!     set velocity grid ocean mask
!-----------------------------------------------------------------------

      umsk(:,:) = 0.
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)
!     remove isolated bays
      do j=2,jmtm1
        do i=2,imtm1
          tmsk(i,j) = max (umsk(i,j), umsk(i-1,j), umsk(i,j-1)
     &,                    umsk(i-1,j-1))
        enddo
      enddo
      call embmbc (tmsk)
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (tmsk(i,j), tmsk(i+1,j), tmsk(i,j+1)
     &,                    tmsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)

      do j=1,jmt
        do i=1,imt
          if (tmsk(i,j) .ge. 0.5) then
            if (hice(i,j,1) .le. 0.) aice(i,j,1) = 0.
            if (hice(i,j,2) .le. 0.) aice(i,j,2) = 0.
          endif
        enddo
      enddo
      call embmbc (aice(1,1,1))
      call embmbc (aice(1,1,2))

!-----------------------------------------------------------------------
!     set the river model
!-----------------------------------------------------------------------

      call rivinit

!-----------------------------------------------------------------------
!     set ocean coalbedo
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .gt. 0) then
!           varies from 0.895 at the equator to 0.815 at the pole
            sbc(i,j,isca) = 0.87 + 0.02*cos(abs(tlat(i,j))*2./radian)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     read vegetation class
!-----------------------------------------------------------------------

      rveg(:,:) = 0.
      fname = new_file_name ("L_potveg.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        call getvara ('L_potveg', iou, imtm2*jmtm2, ib, ic, tmpij
     &,   c1, c0)
        rveg(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (rveg)
      endif
      do j=1,jmt
        do i=1,imt
          iveg(i,j) = iice
          if (rveg(i,j) .gt. 0.6 .and. rveg(i,j) .lt. 7.4)
     &      iveg(i,j) = nint(rveg(i,j))
        enddo
      enddo

      call gvsbc

!----------------------------------------------------------------------
!     initialize elastic viscous plastic variables
!-----------------------------------------------------------------------

      dlam = dxu(int(imt/2))/100.
      dphi = dyu(int(jmt/2))/100.
      diff1 = 0.004
      diff1 = diff1*dlam
      diff2 = diff1*dlam**2
      eccice = 2.
      ecc2 = 1./(eccice**2)
      ecc2m = 2.*(1.-ecc2)
      ecc2p = (1.+ecc2)
      zetamin = 4.e11
      eyc = 0.25
      dte = dtatm/float(ndte)
      dtei = 1./dte
      floor = 1.e-11
      do j=2,jmtm1
        do i=2,imtm1
           xyminevp = (min(cst(j)*dxt(i),dyt(j)))**2
        enddo
      enddo

!-----------------------------------------------------------------------
!     check ice velocity calculation
!-----------------------------------------------------------------------

      if (nivts .gt. nint(segtim*daylen/dtatm)) then
        write(*,*) '==> Warning: ice velocities will be calculated'
        write(*,*) '             every coupling time.'
        nivts =  nint(segtim*daylen/dtatm)
      endif
!-----------------------------------------------------------------------
!     zero time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_tavg (is, ie, js, je, 0)

!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_tsi (is, ie, js, je, 0)

      return
      end
