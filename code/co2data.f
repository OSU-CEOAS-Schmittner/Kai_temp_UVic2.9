! source file: /data/home/kai/dev/UVic2.9/updates/co2data.F
      subroutine co2emitdata

      return
      end

      subroutine co2ccndata
!=======================================================================
!     routine to read and interpolate one dimensional forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer iou, n, ln, ib(10), ic(10)

      logical inqvardef, exists, track

      real avg_co2, dat(3), data_time, pk, tim(3), wt1, wt3, fa

      real, allocatable :: data(:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "cembm.h"
      include "switch.h"
      include "tmngr.h"
      include "atm.h"

!     fa is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      fa = 1./(4.138e-7*rhoatm*shc)
      track = .false.

      if (.not. allocated (time)) then
        name = "A_co2.nc"
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Error: ", trim(fname), " does not exist."
          stop '=>co2ccn'
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          text = 'years'
          call getatttext (iou, 'time', 'units', text)
          if (trim(text) .eq. "days since 1-1-1")
     &      time(:) = time(:)/yrlen - 1.
          if (trim(text) .eq. "days since 0-1-1")
     &       time(:) = time(:)/yrlen
          if (trim(text) .eq. "years since 1-1-1")
     &      time(:) = time(:) - 1.
          exists = inqvardef('A_co2', iou)
          if (.not. exists) then
            print*, "==>  Warning: A_co2 data does not exist."
          else
            call getvara ('A_co2', iou, ln, ib, ic, data, c1, c0)
          endif
         endif
        tim(:) = time(1)
        dat(:) = data(1)
      endif

      tim(2) = min(time(ln), max(time(1), co2_yr))

      if (tim(2) .le. time(1)) then
        dat(2) = data(1)
      elseif (tim(2) .ge. time(ln)) then
        dat(2) = data(ln)
      else
        if (tim(2) .gt. tim(3)) then
          do n=2,ln
            if (time(n-1) .le. tim(2) .and. time(n) .ge. tim(2)) then
              tim(1) = time(n-1)
              dat(1) = data(n-1)
              tim(3) = time(n)
              dat(3) = data(n)
            endif
          enddo
        endif
        wt1 = 1.
        if (tim(3) .ne. tim(1)) wt1 = (tim(3)-tim(2))/(tim(3)-tim(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        dat(2) = dat(1)*wt1 + dat(3)*wt3
      endif

      co2ccn = dat(2)

      return
      end

      subroutine satdata

      return
      end
