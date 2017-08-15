! source file: /data/home/kai/dev/UVic2.9/updates/embm_rest.F
      subroutine embm_rest_in (fname, ids, ide, jds, jde)

!=======================================================================
!     input routine for atmospheric restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================

      implicit none

      character(*) :: fname
      character(32) :: nstamp
      character(3) :: a3
      character(120) :: var1, var2

      integer iou, ln, n, ntrec, ids, ide, jds, jde, ig
      integer ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical inqvardef, exists, use_sealev_data, use_ice_data

      real tmp
      real, allocatable :: tmpij(:,:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "cembm.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "ice.h"
      include "evp.h"
      include "tmngr.h"
      include "switch.h"

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global read domain)
!-----------------------------------------------------------------------
      ils = max(ids,1)
      ile = min(ide,imt)
      jls = max(jds,1)
      jle = min(jde,jmt)

      allocate ( tmpij(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     read 1d data (t)
!-----------------------------------------------------------------------
      tmp = nats
      call getvars ('nats', iou, ntrec, tmp, c1, c0)
      nats = tmp
      tmp = dayoyr
      call getvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      dayoyr = tmp
      tmp = itt
      call getvars ('itt', iou, ntrec, tmp, c1, c0)
      itt = tmp
      tmp = irstdy
      call getvars ('irstdy', iou, ntrec, tmp, c1, c0)
      irstdy = tmp
      tmp = msrsdy
      call getvars ('msrsdy', iou, ntrec, tmp, c1, c0)
      msrsdy = tmp
      tmp = relyr
      call getvars ('relyr', iou, 1, tmp, c1, c0)
      relyr = tmp
      tmp = totaltime
      call getvars ('totaltime', iou, ntrec, tmp, c1, c0)
      totaltime = tmp
      tmp = year0
      call getvars ('year', iou, ntrec, tmp, c1, c0)
      nyear = tmp
      tmp = month0
      call getvars ('month', iou, ntrec, tmp, c1, c0)
      nmonth = tmp
      tmp = day0
      call getvars ('day', iou, ntrec, tmp, c1, c0)
      nday = tmp
      tmp = hour0
      call getvars ('hour', iou, ntrec, tmp, c1, c0)
      nhour = tmp
      tmp = min0
      call getvars ('minute', iou, ntrec, tmp, c1, c0)
      nmin = tmp
      tmp = sec0
      call getvars ('second', iou, ntrec, tmp, c1, c0)
      nsec = tmp
      call mkstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = carbemit
      call getvars ('carbemit', iou, 1, tmp, c1, c0)
      carbemit = tmp
      if (init_time_in) then
        itt = 0
        irstdy = 0
        msrsdy = 0
        relyr = 0.0
        call mkstmp (stamp, year0, month0, day0, hour0, min0, sec0)
        carbemit = 0.
      endif
      use_ice_data = .true.
      use_sealev_data = .true.

!-----------------------------------------------------------------------
!     read 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        var1 = 'at1_'//trim(a3)
        var2 = 'at2_'//trim(a3)
        if (trim(mapat(n)) .eq. 'sat') then
          if (inqvardef('slat1', iou)) var1 = 'slat1'
          if (inqvardef('slat2', iou)) var2 = 'slat2'
        elseif (trim(mapat(n)) .eq. 'shum') then
          if (inqvardef('shum1', iou)) var1 = 'shum1'
          if (inqvardef('shum2', iou)) var2 = 'shum2'
        elseif (trim(mapat(n)) .eq. 'co2') then
          if (inqvardef('co21', iou)) var1 = 'co21'
          if (inqvardef('co22', iou)) var2 = 'co22'
        endif
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,1,n)
        call getvara(trim(var1), iou, ln, ib, ic, tmpij, c1, c0)
        at(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,2,n)
        call getvara(trim(var2), iou, ln, ib, ic, tmpij, c1, c0)
        at(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
      enddo
      tmpij(ils:ile,jls:jle) = rh(ils:ile,jls:jle)
      call getvara ('rh', iou, ln, ib, ic, tmpij, c1, c0)
      rh(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call getvara ('precip', iou, ln, ib, ic, tmpij, c1, c0)
      precip(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isst)
      call getvara ('sbc_sst', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isst) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isss)
      call getvara ('sbc_sss', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isss) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ihflx)
      call getvara ('sbc_hflx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ihflx) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isflx)
      call getvara ('sbc_sflx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isflx) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic)
      call getvara ('sbc_ssdic', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdic) = tmpij(ils:ile,jls:jle)
       tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic13)
      call getvara ('sbc_ssdic13', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdic13) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isspo4)
      call getvara ('sbc_sspo4', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isspo4) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isso2)
      call getvara ('sbc_sso2', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isso2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issalk)
      call getvara ('sbc_ssalk', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issalk) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdop)
      call getvara ('sbc_ssdop', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdop) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issno3)
      call getvara ('sbc_ssno3', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issno3) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdon)
      call getvara ('sbc_ssdon', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdon) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdin15)
      call getvara ('sbc_ssdin15', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdin15) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdon15)
      call getvara ('sbc_ssdon15', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdon15) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdfe)
      call getvara ('sbc_ssdfe', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdfe) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdoc13)
      call getvara ('sbc_ssdoc13', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdoc13) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call getvara ('rtbar', iou, ln, ib, ic, tmpij, c1, c0)
      rtbar(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = atbar(ils:ile,jls:jle)
      call getvara ('atbar', iou, ln, ib, ic, tmpij, c1, c0)
      atbar(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,1)
      call getvara ('soilm1', iou, ln, ib, ic, tmpij, c1, c0)
      soilm(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,2)
      call getvara ('soilm2', iou, ln, ib, ic, tmpij, c1, c0)
      soilm(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call getvara ('surf', iou, ln, ib, ic, tmpij, c1, c0)
      surf(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = tice(ils:ile,jls:jle)
      call getvara ('tice', iou, ln, ib, ic, tmpij, c1, c0)
      tice(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,1)
      call getvara ('hice1', iou, ln, ib, ic, tmpij, c1, c0)
      hice(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,2)
      call getvara ('hice2', iou, ln, ib, ic, tmpij, c1, c0)
      hice(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle,1)
      call getvara ('aice1', iou, ln, ib, ic, tmpij, c1, c0)
      aice(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle,2)
      call getvara ('aice2', iou, ln, ib, ic, tmpij, c1, c0)
      aice(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle,1)
      call getvara ('hsno1', iou, ln, ib, ic, tmpij, c1, c0)
      hsno(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle,2)
      call getvara ('hsno2', iou, ln, ib, ic, tmpij, c1, c0)
      hsno(ils:ile,jls:jle,2)= tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call getvara ('uice', iou, ln, ib, ic, tmpij, c1, c0)
      uice(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call getvara ('vice', iou, ln, ib, ic, tmpij, c1, c0)
      vice(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isu)
      call getvara ('su', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isu) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isv)
      call getvara ('sv', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isv) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igu)
      call getvara ('gu', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,igu) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igv)
      call getvara ('gv', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,igv) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11n(ils:ile,jls:jle)
      call getvara ('sig11n', iou, ln, ib, ic, tmpij, c1, c0)
      sig11n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11s(ils:ile,jls:jle)
      call getvara ('sig11s', iou, ln, ib, ic, tmpij, c1, c0)
      sig11s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11e(ils:ile,jls:jle)
      call getvara ('sig11e', iou, ln, ib, ic, tmpij, c1, c0)
      sig11e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11w(ils:ile,jls:jle)
      call getvara ('sig11w', iou, ln, ib, ic, tmpij, c1, c0)
      sig11w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22n(ils:ile,jls:jle)
      call getvara ('sig22n', iou, ln, ib, ic, tmpij, c1, c0)
      sig22n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22s(ils:ile,jls:jle)
      call getvara ('sig22s', iou, ln, ib, ic, tmpij, c1, c0)
      sig22s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22e(ils:ile,jls:jle)
      call getvara ('sig22e', iou, ln, ib, ic, tmpij, c1, c0)
      sig22e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22w(ils:ile,jls:jle)
      call getvara ('sig22w', iou, ln, ib, ic, tmpij, c1, c0)
      sig22w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12n(ils:ile,jls:jle)
      call getvara ('sig12n', iou, ln, ib, ic, tmpij, c1, c0)
      sig12n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12s(ils:ile,jls:jle)
      call getvara ('sig12s', iou, ln, ib, ic, tmpij, c1, c0)
      sig12s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12e(ils:ile,jls:jle)
      call getvara ('sig12e', iou, ln, ib, ic, tmpij, c1, c0)
      sig12e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12w(ils:ile,jls:jle)
      call getvara ('sig12w', iou, ln, ib, ic, tmpij, c1, c0)
      sig12w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)

      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Atm restart read from ',trim(fname),' on ', nstamp

      deallocate ( tmpij )

      return
      end

      subroutine embm_rest_def (fname)
!=======================================================================
!     definition routine for atmospheric restarts

!   inputs:
!     fname = file name
!=======================================================================

      implicit none

      character(*) :: fname
      character(3) :: a3

      integer iou, n, igs, ige, ig, jgs, jge, jg, it(10), iu(10)
      integer id_time, id_xt, id_xu, id_yt, id_yu, id_xt_e, id_xu_e
      integer id_yt_e, id_yu_e, id_track

      real c1e20

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "ice.h"
      include "evp.h"
      include "iounit.h"
      include "tmngr.h"

      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('longitude', iou, ig, id_xt)
      call defdim ('latitude', iou, jg, id_yt)
      call defdim ('longitude_V', iou, ig, id_xu)
      call defdim ('latitude_V', iou, jg, id_yu)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('longitude_V_edges', iou, ig+1, id_xu_e)
      call defdim ('latitude_V_edges', iou, jg+1, id_yu_e)
      call defdim ('track', iou, mtrack, id_track)

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'years since 0-1-1')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('nats', iou, 1, it, c0, c0, ' ', 'D'
     &, 'nats', ' ',' ')
      call defvar ('dayoyr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'dayoyr', ' ',' ')
      call defvar ('itt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itt', ' ',' ')
      call defvar ('irstdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'irstdy', ' ',' ')
      call defvar ('msrsdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'msrsdy', ' ',' ')
      call defvar ('relyr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'relyr', ' ',' ')
      call defvar ('totaltime', iou, 1, it, c0, c0, ' ', 'D'
     &, 'totaltime', ' ',' ')
      call defvar ('year', iou, 1, it, c0, c0, ' ', 'D'
     &, 'year', ' ',' ')
      call defvar ('month', iou, 1, it, c0, c0, ' ', 'D'
     &, 'month', ' ',' ')
      call defvar ('day', iou, 1, it, c0, c0, ' ', 'D'
     &, 'day', ' ',' ')
      call defvar ('hour', iou, 1, it, c0, c0, ' ', 'D'
     &, 'hour', ' ',' ')
      call defvar ('minute', iou, 1, it, c0, c0, ' ', 'D'
     &, 'minute', ' ',' ')
      call defvar ('second', iou, 1, it, c0, c0, ' ', 'D'
     &, 'second', ' ',' ')
      call defvar ('carbemit', iou, 1, it, c0, c0, ' ', 'D'
     &, 'carbemit', ' ',' ')
      call defvar ('co2ccn', iou, 1, it, c0, c0, ' ', 'D'
     &, 'co2ccn', ' ',' ')
      call defvar ('c13ccn', iou, 1, it, c0, c0, ' ', 'D'
     &, 'c13ccn', ' ',' ')
      call defvar ('ice_yr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'ice_yr', ' ',' ')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yt
      call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_xu
      call defvar ('longitude_V', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yu
      call defvar ('latitude_V', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('slat1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'sea level atmospheric temperature at tau', ' ', ' ')
          call defvar ('slat2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'sea level atmospheric temperature at tau+1', ' ', ' ')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('shum1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric surface specific humidity at tau', ' ', ' ')
          call defvar ('shum2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric surface specific humidity at tau+1', ' ', ' ')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('co21', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric co2 at tau', ' ', ' ')
          call defvar ('co22', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric co2 at tau+1', ' ', ' ')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('at1_'//trim(a3), iou , 3, it, -c1e20, c1e20, ' '
     &,     'D', 'at1_'//trim(a3), ' ', ' ')
          call defvar ('at2_'//trim(a3), iou , 3, it, -c1e20, c1e20, ' '
     &,     'D', 'at2_'//trim(a3), ' ', ' ')
        endif
      enddo
      call defvar ('rh', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'rh', ' ', ' ')
      call defvar ('precip', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'precip', ' ', ' ')
      call defvar ('sbc_sst', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sst', ' ', ' ')
      call defvar ('sbc_sss', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sss', ' ', ' ')
      call defvar ('sbc_hflx', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sst', ' ', ' ')
      call defvar ('sbc_sflx', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sss', ' ', ' ')
      call defvar ('sbc_ssdic', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdic', ' ', ' ')
      call defvar ('sbc_ssdic13', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdic13', ' ', ' ')
      call defvar ('sbc_sspo4', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sspo4', ' ', ' ')
      call defvar ('sbc_sso2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sso2', ' ', ' ')
      call defvar ('sbc_ssalk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssalk', ' ', ' ')
      call defvar ('sbc_ssdop', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdop', ' ', ' ')
      call defvar ('sbc_ssno3', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssno3', ' ', ' ')
      call defvar ('sbc_ssdon', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdon', ' ', ' ')
      call defvar ('sbc_ssdin15', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdin15', ' ', ' ')
      call defvar ('sbc_ssdon15', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdon15', ' ', ' ')
      call defvar ('sbc_ssdfe', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdfe', ' ', ' ')
      call defvar ('sbc_ssdoc13', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdoc13', ' ', ' ')
      call defvar ('tmsk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'tmsk', ' ', ' ')
      call defvar ('hicel', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hicel', ' ', ' ')
      call defvar ('aicel', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aicel', ' ', ' ')
      call defvar ('rtbar', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'rtbar', ' ', ' ')
      call defvar ('atbar', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'atbar', ' ', ' ')
      call defvar ('soilm1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'soilm1', ' ', ' ')
      call defvar ('soilm2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'soilm2', ' ', ' ')
      call defvar ('surf', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'surf', ' ', ' ')
      call defvar ('tice', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'tice', ' ', ' ')
      call defvar ('hice1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hice1', ' ', ' ')
      call defvar ('hice2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hice2', ' ', ' ')
      call defvar ('aice1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aice1', ' ', ' ')
      call defvar ('aice2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aice2', ' ', ' ')
      call defvar ('hsno1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsno1', ' ', ' ')
      call defvar ('hsno2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsno2', ' ', ' ')
      call defvar ('uice', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'uice', ' ', ' ')
      call defvar ('vice', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'vice', ' ', ' ')
      call defvar ('su', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'su', ' ', ' ')
      call defvar ('sv', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'sv', ' ', ' ')
      call defvar ('gu', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'gu', ' ', ' ')
      call defvar ('gv', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'gv', ' ', ' ')
      call defvar ('sig11n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11n', ' ', ' ')
      call defvar ('sig11s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11s', ' ', ' ')
      call defvar ('sig11e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11e', ' ', ' ')
      call defvar ('sig11w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11w', ' ', ' ')
      call defvar ('sig22n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22n', ' ', ' ')
      call defvar ('sig22s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22s', ' ', ' ')
      call defvar ('sig22e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22e', ' ', ' ')
      call defvar ('sig22w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22w', ' ', ' ')
      call defvar ('sig12n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12n', ' ', ' ')
      call defvar ('sig12s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12s', ' ', ' ')
      call defvar ('sig12e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12e', ' ', ' ')
      call defvar ('sig12w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12w', ' ', ' ')
      call enddef (iou)

      return
      end

      subroutine embm_rest_out (fname, ids, ide, jds, jde)
!=======================================================================
!     output routine for atmospheric restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================

      implicit none

      character(*) :: fname
      character(3) :: a3
      character(32) :: nstamp
      character(120) :: var1, var2

      integer i, iou, j, ln, n, ntrec, ids, ide, jds, jde, igs, ige, ig
      integer jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical inqvardef, exists

      real time, tmp
      real, allocatable :: tmpij(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "cembm.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "ice.h"
      include "evp.h"
      include "tmngr.h"
      include "iounit.h"
      include "switch.h"
      include "npzd.h"
      real xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1), yu_e(jmt+1)

      nstamp = stamp

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      if (init_time_out) then
        tmp = 0.
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('relyr', iou, ntrec, tmp, c1, c0)
        call mkstmp (nstamp, year0, month0, day0, hour0, min0, sec0)
        carbemit = 0.
      else
        tmp = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = itt
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = iday(imodeltime)
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = msday(imodeltime)
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
        tmp = relyr
        call putvars ('relyr', iou, ntrec, tmp, c1, c0)
      endif
      tmp = nats
      call putvars ('nats', iou, ntrec, tmp, c1, c0)
      tmp = dayoyr
      call putvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      call putvars ('totaltime', iou, ntrec, totaltime, c1, c0)
      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, ntrec, tmp, c1, c0)
      tmp = carbemit
      call putvars ('carbemit', iou, ntrec, tmp, c1, c0)
      tmp = co2ccn
      call putvars ('co2ccn', iou, ntrec, tmp, c1, c0)
      c13ccn = (1 + dc13ccn*0.001)*rc13std*co2ccn
     &      /(1+(1 + dc13ccn*0.001)*rc13std)
      tmp = c13ccn
      print*,'write c13ccn=',c13ccn
      call putvars ('c13ccn', iou, ntrec, tmp, c1, c0)
      tmp = ice_yr
      call putvars ('ice_yr', iou, ntrec, tmp, c1, c0)

!-----------------------------------------------------------------------
!     write 1d data (x or y)
!-----------------------------------------------------------------------
      allocate ( tmpi(igs:ige) )
      allocate ( tmpj(jgs:jge) )
      allocate ( tmpie(igs:ige+1) )
      allocate ( tmpje(jgs:jge+1) )

      ib(1) = 1
      ic(1) = ig
      tmpi(igs:ige) = xt(igs:ige)
      call putvara ('longitude', iou, ig, ib, ic, tmpi, c1, c0)
      tmpi(igs:ige) = xu(igs:ige)
      call putvara ('longitude_V', iou, ig, ib, ic, tmpi, c1, c0)

      ic(1) = jg
      tmpj(jgs:jge) = yt(jgs:jge)
      call putvara ('latitude', iou, jg, ib, ic, tmpj, c1, c0)
      tmpj(jgs:jge) = yu(jgs:jge)
      call putvara ('latitude_V', iou, jg, ib, ic, tmpj, c1, c0)

      ic(1) = ig + 1
      call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xt_e(igs:ige+1)
      call putvara ('longitude_edges', iou, ig+1, ib, ic, tmpie
     &, c1, c0)
      call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xu_e(igs:ige+1)
      call putvara ('longitude_V_edges', iou, ig+1, ib, ic, tmpie
     &, c1, c0)

      ic(1) = jg + 1
      call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
      call putvara ('latitude_edges', iou, jg+1, ib, ic, tmpje
     &, c1, c0)
      call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
      call putvara ('latitude_V_edges', iou, jg+1, ib, ic, tmpje
     &, c1, c0)

      deallocate ( tmpi )
      deallocate ( tmpj )
      deallocate ( tmpie )
      deallocate ( tmpje )

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          var1 = 'slat1'
          var2 = 'slat2'
        elseif (trim(mapat(n)) .eq. 'shum') then
          var1 = 'shum1'
          var2 = 'shum2'
        elseif (trim(mapat(n)) .eq. 'co2') then
          var1 = 'co21'
          var2 = 'co22'
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          var1 = 'at1_'//trim(a3)
          var2 = 'at2_'//trim(a3)
        endif
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,1,n)
        call putvara(trim(var1), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,2,n)
        call putvara(trim(var2), iou, ln, ib, ic, tmpij, c1, c0)
      enddo
      tmpij(ils:ile,jls:jle) = rh(ils:ile,jls:jle)
      call putvara ('rh', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call putvara ('precip', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isst)
      call putvara ('sbc_sst', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isss)
      call putvara ('sbc_sss', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ihflx)
      call putvara ('sbc_hflx', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isflx)
      call putvara ('sbc_sflx', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic)
      call putvara ('sbc_ssdic', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic13)
      call putvara ('sbc_ssdic13', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isspo4)
      call putvara ('sbc_sspo4', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isso2)
      call putvara ('sbc_sso2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issalk)
      call putvara ('sbc_ssalk', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdop)
      call putvara ('sbc_ssdop', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issno3)
      call putvara ('sbc_ssno3', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdon)
      call putvara ('sbc_ssdon', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdin15)
      call putvara ('sbc_ssdin15', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdon15)
      call putvara ('sbc_ssdon15', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdfe)
      call putvara ('sbc_ssdfe', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdoc13)
      call putvara ('sbc_ssdoc13', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) =  tmsk(ils:ile,jls:jle)
      call putvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hicel(ils:ile,jls:jle,2)
      call putvara ('hicel', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = aicel(ils:ile,jls:jle,2)
      call putvara ('aicel', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call putvara ('rtbar', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = atbar(ils:ile,jls:jle)
      call putvara ('atbar', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,1)
      call putvara ('soilm1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,2)
      call putvara ('soilm2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call putvara ('surf', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = tice(ils:ile,jls:jle)
      call putvara ('tice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,1)
      call putvara ('hice1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,2)
      call putvara ('hice2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle,1)
      call putvara ('aice1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle,2)
      call putvara ('aice2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle,1)
      call putvara ('hsno1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle,2)
      call putvara ('hsno2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call putvara ('uice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call putvara ('vice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isu)
      call putvara ('su', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isv)
      call putvara ('sv', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igu)
      call putvara ('gu', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igv)
      call putvara ('gv', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11n(ils:ile,jls:jle)
      call putvara ('sig11n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11s(ils:ile,jls:jle)
      call putvara ('sig11s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11e(ils:ile,jls:jle)
      call putvara ('sig11e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11w(ils:ile,jls:jle)
      call putvara ('sig11w', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22n(ils:ile,jls:jle)
      call putvara ('sig22n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22s(ils:ile,jls:jle)
      call putvara ('sig22s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22e(ils:ile,jls:jle)
      call putvara ('sig22e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22w(ils:ile,jls:jle)
      call putvara ('sig22w', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12n(ils:ile,jls:jle)
      call putvara ('sig12n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12s(ils:ile,jls:jle)
      call putvara ('sig12s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12e(ils:ile,jls:jle)
      call putvara ('sig12e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12w(ils:ile,jls:jle)
      call putvara ('sig12w', iou, ln, ib, ic, tmpij, c1, c0)

      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Atm restart written to ',trim(fname),' on ', nstamp

      deallocate ( tmpij )

      return
      end
