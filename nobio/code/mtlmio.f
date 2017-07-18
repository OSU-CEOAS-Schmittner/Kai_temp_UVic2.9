! source file: /data/home/kai/dev/UVic2.9/nobio/updates/mtlmio.F
      subroutine mtlmout (is, ie, js, je)

!-----------------------------------------------------------------------
!     Output routine for the mtlm
!-----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "coord.h"
      include "grdvar.h"
      include "mtlm.h"
      include "csbc.h"
      include "cembm.h"
      include "iounit.h"
      include "switch.h"
      include "tmngr.h"

      character(120) :: fname
      character(32) :: nstamp

      integer is, ie, js, je, ntrec
      integer nyear, nmonth, nday, nhour, nmin, nsec

      real avgper, time, tmp

      if (tsits .and. ntatil .ne. 0) then

        call ta_mtlm_tsi (is, ie, js, je, 2)

        avgper = tsiper*accel
        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        call def_tsi
        call def_tsi_mtlm (fname)
        tai_clnd = 0.
        tai_cfa2l = 0.
        avgper = tsiper*accel
        if (avgper .le. 1e-6) avgper = 0.
        tmp = 0.5
        time = time - tmp*avgper/365.

        call mtlm_tsi_out (fname, avgper, time, nstamp, tai_CS
     &,   tai_RESP_S, tai_LIT_C_T, tai_BURN, tai_CV, tai_NPP, tai_GPP
     &,   tai_HT, tai_LAI, tai_LYING_SNOW, tai_TSOIL, tai_TSTAR
     &,   tai_M, tai_ET, tai_clnd, tai_cfa2l, ntrec
     &    )
        call ta_mtlm_tsi (is, ie, js, je, 0)

      endif
      if (timavgts .and. ntatsl .ne. 0) then

!-----------------------------------------------------------------------
!       write atmospheric time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_mtlm_tavg (is, ie, js, je, 2)

!       write time averaged data

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        call def_tavg
        call def_tavg_mtlm (fname)
        avgper = timavgper*accel
        if (avgper .le. 1e-6) avgper = 0.
        tmp = 0.5
        time = time - tmp*avgper/365.

        call mtlm_tavg_out (fname, is, ie, js, je, imt, jmt
     &,   POINTS, NPFT, NTYPE, xt, yt, xu, yu, dxt, dyt, dxu, dyu
     &,   avgper, time, nstamp, land_map, ta_TS1, ta_CS, ta_RESP_S
     &,   ta_LIT_C_T, ta_BURN, ta_FRAC, ta_GPP, ta_NPP, ta_HT, ta_LAI
     &,   ta_C_VEG
     &,   tlat, tlon, tgarea, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Lnd time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

!       zero time average accumulators

        call ta_mtlm_tavg (is, ie, js, je, 0)

      endif

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_mtlm (0, fname)
          call mtlm_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_mtlm (1, fname)
          call mtlm_rest_out (fname, is, ie, js, je)
        endif
      endif

      return
      end

      subroutine ta_mtlm_tavg (is, ie, js, je, iflag)

!=======================================================================
!     land data time averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      include "size.h"
      include "mtlm.h"
      include "csbc.h"
      include "switch.h"

      integer i, is, ie, j, js, je, iflag, L, n

      real rntatsl

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatsl = 0
        ta_TS1(:) = 0.
        ta_TSTAR_GB(:) = 0.
        ta_M(:) = 0.
        ta_CS(:) = 0.
        ta_RESP_S(:) = 0.
        ta_LIT_C_T(:) = 0.
        ta_BURN(:) = 0.
        ta_GPP(:,:) = 0.
        ta_NPP(:,:) = 0.
        ta_HT(:,:) = 0.
        ta_LAI(:,:) = 0.
        ta_C_VEG(:,:) = 0.
        ta_LYING_SNOW(:) = 0.
        ta_SURF_ROFF(:) = 0.
        ta_FRAC(:,:) = 0.

      elseif (iflag .eq. 1) then

!       accumulate
        ntatsl = ntatsl + 1
        ta_TS1(:) = ta_TS1(:) + TS1(:)
        ta_TSTAR_GB(:) = ta_TSTAR_GB(:) + TSTAR_GB(:)
        ta_M(:) = ta_M(:) + M(:)
        ta_CS(:) = ta_CS(:) + CS(:)
        ta_RESP_S(:) = ta_RESP_S(:) + RESP_S(:)
        ta_LIT_C_T(:) = ta_LIT_C_T(:) + LIT_C_T(:)

        do n=1,NPFT
          ta_GPP(:,n) = ta_GPP(:,n) + GPP(:,n)*FRAC(:,n)
          ta_NPP(:,n) = ta_NPP(:,n) + NPP(:,n)*FRAC(:,n)
          ta_HT(:,n) = ta_HT(:,n) + HT(:,n)*FRAC(:,n)
          ta_LAI(:,n) = ta_LAI(:,n) + LAI(:,n)*FRAC(:,n)
          ta_C_VEG(:,n) = ta_C_VEG(:,n) + C_VEG(:,n)*FRAC(:,n)
        enddo

        ta_LYING_SNOW(:) = ta_LYING_SNOW(:) + LYING_SNOW(:)
        ta_SURF_ROFF(:) = ta_SURF_ROFF(:) + SURF_ROFF(:)
        ta_FRAC(:,:) = ta_FRAC(:,:) + FRAC(:,:)

      elseif (iflag .eq. 2 .and. ntatsl .ne. 0) then

!       average
        rntatsl = 1./float(ntatsl)
        ta_TS1(:) = ta_TS1(:)*rntatsl
        ta_TSTAR_GB(:) = ta_TSTAR_GB(:) *rntatsl
        ta_M(:) = ta_M(:)*rntatsl
        ta_CS(:) = ta_CS(:)*rntatsl
        ta_RESP_S(:) = ta_RESP_S(:)*rntatsl
        ta_LIT_C_T(:) = ta_LIT_C_T(:)*rntatsl/SEC_YEAR
        ta_BURN(:) = ta_BURN(:)*rntatsl
        ta_GPP(:,:) = ta_GPP(:,:)*rntatsl
        ta_NPP(:,:) = ta_NPP(:,:)*rntatsl
        ta_HT(:,:) = ta_HT(:,:)*rntatsl
        ta_LAI(:,:) = ta_LAI(:,:)*rntatsl
        ta_C_VEG(:,:) = ta_C_VEG(:,:)*rntatsl
        ta_LYING_SNOW(:) = ta_LYING_SNOW(:)*rntatsl
        ta_SURF_ROFF(:) = ta_SURF_ROFF(:)*rntatsl
        ta_FRAC(:,:) = ta_FRAC(:,:)*rntatsl
      endif

      return
      end

      subroutine ta_mtlm_tsi (is, ie, js, je, iflag)

!=======================================================================
!     land data time integral averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      include "size.h"
      include "csbc.h"
      include "mtlm.h"
      include "switch.h"

      integer is, ie, js, je, iflag, n

      real rntatil, data(imt,jmt), dmsk(imt,jmt), wt(imt, jmt), tmp

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

       if (iflag .eq. 0.) then

!       zero
        ntatil = 0
        tai_CS = 0
        tai_RESP_S = 0
        tai_LIT_C_T = 0
        tai_BURN = 0
        tai_CV = 0
        tai_NPP = 0
        tai_GPP = 0
        tai_HT = 0
        tai_LAI = 0
        tai_LYING_SNOW = 0
        tai_TSOIL = 0
        tai_TSTAR = 0
        tai_M = 0
        tai_ET = 0
      elseif (iflag .eq. 1) then

!       set data mask
        dmsk(:,:) = 1.
        where (land_map(:,:) .eq. 0) dmsk(:,:) = 0.
!       accumulate
        ntatil = ntatil + 1

        call unloadland (POINTS, CS, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_CS = tai_CS + tmp*1.e-4

        call unloadland (POINTS, RESP_S, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_RESP_S = tai_RESP_S + tmp*1.e-4

        call unloadland (POINTS, LIT_C_T, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_LIT_C_T = tai_LIT_C_T + tmp*1.e-4

        call unloadland (POINTS, CV, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_CV = tai_CV + tmp*1.e-4

        do n=1,NPFT

          call unloadland (POINTS, FRAC(1,n), imt, jmt, land_map, wt)

          call unloadland (POINTS, NPP(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areatot (data, dmsk, tmp)
!         convert area to m2
          tai_NPP = tai_NPP + tmp*1.e-4

          call unloadland (POINTS, GPP(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areatot (data, dmsk, tmp)
!         convert area to m2
          tai_GPP = tai_GPP + tmp*1.e-4

          call unloadland (POINTS, HT(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areaavg (data, dmsk, tmp)
          tai_HT = tai_HT + tmp

          call unloadland (POINTS, LAI(1,n), imt, jmt, land_map, data)
          data(:,:) = data(:,:)*wt(:,:)
          call areaavg (data, dmsk, tmp)
          tai_LAI = tai_LAI + tmp

        enddo

        call unloadland (POINTS, LYING_SNOW, imt, jmt, land_map, data)
        call areatot (data, dmsk, tmp)
!       convert area to m2
        tai_LYING_SNOW = tai_LYING_SNOW + tmp*1.e-4

        call unloadland (POINTS, TS1, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_TSOIL = tai_TSOIL + tmp

        call unloadland (POINTS, TSTAR_GB, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_TSTAR = tai_TSTAR + tmp

        call unloadland (POINTS, M, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_M = tai_M + tmp

        call unloadland (POINTS, ET, imt, jmt, land_map, data)
        call areaavg (data, dmsk, tmp)
        tai_ET = tai_ET + tmp

      elseif (iflag .eq. 2 .and. ntatil .ne. 0) then

!       average
        rntatil = 1./float(ntatil)
        tai_CS = tai_CS*rntatil
        tai_RESP_S = tai_RESP_S*rntatil
        tai_LIT_C_T = tai_LIT_C_T*rntatil/SEC_YEAR
        tai_BURN = tai_BURN*rntatil
        tai_CV = tai_CV*rntatil
        tai_NPP = tai_NPP*rntatil
        tai_GPP = tai_GPP*rntatil
        tai_HT = tai_HT*rntatil
        tai_LAI = tai_LAI*rntatil
        tai_LYING_SNOW = tai_LYING_SNOW*rntatil
        tai_TSOIL = tai_TSOIL*rntatil
        tai_TSTAR = tai_TSTAR*rntatil
        tai_M = tai_M*rntatil
        tai_ET = tai_ET*rntatil
      endif

      return
      end

      subroutine unloadland (ld, dl, id, jd, map, dij)

      implicit none

      integer i, id, j, jd, l, ld, map(id,jd)
      real dl(ld), dij(id,jd)

      dij(:,:) = 0.
      do j=1,jd
        do i=1,id
          l = map(i,j)
          if (l .ge. 1 .and. l .le. ld) dij(i,j) = dl(l)
        enddo
      enddo

      return
      end

      subroutine loadland (ld, dl, id, jd, map, dij)

      implicit none

      integer i, id, j, jd, l, ld, map(id,jd)
      real dl(ld), dij(id,jd)

      dl(:) = 0.
      do j=1,jd
        do i=1,id
          l = map(i,j)
          if (l .ge. 1 .and. l .le. ld) dl(l) = dij(i,j)
        enddo
      enddo

      return
      end
