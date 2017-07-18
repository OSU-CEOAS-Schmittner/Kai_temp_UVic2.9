! source file: /data/home/kai/dev/UVic2.9/nobio/updates/def_files.F
      subroutine def_tavg
!=======================================================================
!     defines tavg files for UVic_ESCM
!=======================================================================

      implicit none

      character(120) :: fname

      call def_tavg_embm (fname)
      call def_tavg_mtlm (fname)
      call def_tavg_mom (fname)

      return
      end

      subroutine def_rest (last)
!=======================================================================
!     defines rest files for UVic_ESCM

!     input:
!       last = last restart flag (1 = last)
!=======================================================================

      implicit none

      character(120) :: fname
      integer last

      call def_rest_embm (last, fname)
      call def_rest_mtlm (last, fname)
      call def_rest_mom (last, fname)

      return
      end

      subroutine def_tsi
!=======================================================================
!     defines tsi files for UVic_ESCM
!=======================================================================

      implicit none

      character(120) :: fname

      call def_tsi_embm (fname)
      call def_tsi_mtlm (fname)
      call def_tsi_mom (fname)

      return
      end

      subroutine inqdefined (name, defined)

!=======================================================================
!     keeps track of which files have been defined. assumes the file is
!     about to be defined.

!     input:
!       name = file name

!     output:
!       defined = logical flag
!=======================================================================

      integer n, max_num_files
      parameter (max_num_files=201)

      character(120) :: name
      character(120), allocatable :: file_names(:)

      logical defined

      save file_names

      if (.not. allocated (file_names)) then
        allocate ( file_names(max_num_files) )
        file_names(1:max_num_files) = " "
      endif

      n = 1
      defined = .false.
      do while (.not. defined .and. file_names(n) .ne. " ")
        if (trim(name) .eq. trim(file_names(n))) defined = .true.
        n = n + 1
      enddo
!     always leave file_names(max_num_files) = " " to avoid more
!     testing in the "do while" loop
      if (.not. defined) then
        if (n .lt. max_num_files) then
          file_names(n) = name
        else
          stop 'maximum number of files exceeded in inqdefined'
        endif
      endif

      return
      end

      subroutine def_tavg_embm (fname)
!=======================================================================
!     defines tavg files for the embm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "atm.h"
      include "ice.h"
      include "iounit.h"
      include "tmngr.h"

      defined = .false.
      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_embm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call embm_tavg_def (name, imt, jmt, nat, ncat, xt, yt
     &,                     calendar, expnam, runstamp, mapat)
      endif
      fname = name

      return
      end

      subroutine def_rest_embm (last, fname)
!=======================================================================
!     defines restart file for the embm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_embm.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_embm',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call embm_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_embm (fname)
!=======================================================================
!     defines tsi files for the embm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_embm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call embm_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end

      subroutine def_tavg_mtlm (fname)
!=======================================================================
!     defines tavg files for the mtlm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_mtlm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mtlm_tavg_def (name, imt, jmt, NPFT, NTYPE, xt, yt
     &,                     calendar, expnam, runstamp)
      endif
      fname = name

      return
      end

      subroutine def_rest_mtlm (last, fname)
!=======================================================================
!     defines restart file for the mtlm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_mtlm.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_mtlm',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call mtlm_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_mtlm (fname)
!=======================================================================
!     defines tsi files for the mtlm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_mtlm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mtlm_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end

      subroutine def_tavg_mom (fname)
!=======================================================================
!     defines tavg files for the mom

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "iounit.h"
      include "mw.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_mom',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mom_tavg_def (name, imt, jmt, km, nt, kpzd, xt, yt
     &,                    calendar, expnam, runstamp, mapt)
      endif
      fname = name

      return
      end

      subroutine def_rest_mom (last, fname)
!=======================================================================
!     defines restart file for the mom

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

       implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_mom.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_mom',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call mom_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_mom (fname)
!=======================================================================
!     defines tsi files for the mom

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_mom',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mom_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end

