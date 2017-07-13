! source file: /raid24/aho/UVic2.9/default_comb2/nobio/updates/data.F
      subroutine data (is, ie, js, je)

!----------------------------------------------------------------------
!     data routine
!----------------------------------------------------------------------

      implicit none

      integer ie, is, je, js, i, iou, j, m, n

      real damp1, damp2, realdays, wnext
      real c10, c100, p001, p035

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "ctdbc.h"
      include "tmngr.h"
      include "switch.h"
      include "scalar.h"
      include "atm.h"

      c10 = 10.
      c100 = 100.
      p001 = 0.001
      p035 = 0.035

!-----------------------------------------------------------------------
!     determine the disk pointers, time weight interpolation factor,
!     and whether or not it is time to bring in new S.B.C. from disk
!     based on the time (days) in MOM since dec 31, 1899 midnight.

!     express model time in days after start of S.B.C. by adding time
!     of I.C. to current model time then subtract time at start of
!     S.B.C.. Note that "itemptime" was allocated in settmngr and is
!     only needed as a temporary.
!     need to add "dt" to the model time because the call to
!     data precedes the time stepping loop which calls mom, so the
!     model time has not yet been incremented when data executes.
!-----------------------------------------------------------------------

      do n=1,ntdbc
        call addtime (initial, imodeltime, itemptime)
        call addtime (itemptime, idt, itemptime)
        call subtime (itemptime, isbcstart(n), itemptime)
        daysbc(n) = realdays(itemptime)
      enddo

!-----------------------------------------------------------------------
!     determine the disk pointers, time weight interpolation factor,
!     and whether or not it is time to bring in new S.B.C. from disk
!     based on the time (days) in MOM since dec 31, 1899 midnight.
!-----------------------------------------------------------------------

      do n=1,ntdbc
        call timeinterp (daysbc(n), n, tdrec(1,n), aprec(1,n)
     &,    ntdrec(n), period(n), method, inextd(n), iprevd(n)
     &,    wprev(n), rdtdbc(n), inextm(n), iprevm(n))
        rdtdbc(n) = .false.
        iprevm(n) = iprevd(n)
        inextm(n) = inextd(n)
      enddo

!-----------------------------------------------------------------------
!     read in data for each S.B.C. when necessary
!-----------------------------------------------------------------------

      n = 1
      do m=1,numsbc

        if ( m .eq. itaux ) then
!         x component of windstress
          call get_tdsbc (n, 'O_tau.nc', 'O_tauX', itaux
     &,     rdtdbc(n), c10, c0)

        elseif ( m .eq. itauy ) then
!         y component of windstress
          call get_tdsbc (n, 'O_tau.nc', 'O_tauY', itauy
     &,     rdtdbc(n), c10, c0)

        elseif ( m .eq. iws ) then
!         surface wind speed
          call get_tdsbc (n, 'A_windsur.nc', 'A_windspd', iws
     &,     rdtdbc(n), c100, c0)

        elseif ( m .eq. iaca ) then
!         atmospheric coalbedo
          call get_tdsbc (n, 'A_calbatm.nc', 'A_calbatm', iaca
     &,     rdtdbc(n), c1, c0)

        elseif ( m .eq. iwxq ) then
!         x component of advecting wind
          call get_tdsbc (n, 'A_wind.nc', 'A_windqX', iwxq
     &,     rdtdbc(n), c100, c0)

        elseif ( m .eq. iwyq ) then
!         y component of advecting wind
          call get_tdsbc (n, 'A_wind.nc', 'A_windqY', iwyq
     &,     rdtdbc(n), c100, c0)

        elseif ( m .eq. iwxt ) then
!         x component of advecting wind
          call get_tdsbc (n, 'A_wind.nc', 'A_windtX', iwxt
     &,     rdtdbc(n), c100, c0)

        elseif ( m .eq. iwyt ) then
!         y component of advecting wind
          call get_tdsbc (n, 'A_wind.nc', 'A_windtY', iwyt
     &,     rdtdbc(n), c100, c0)

        elseif ( m .eq. idtr ) then
!         diurnal temperature range
          call get_tdsbc (n, 'L_diurtemp.nc', 'L_diurtemp', idtr
     &,     rdtdbc(n), c1, c0)

        endif

      enddo

      return
      end

      subroutine get_tdsbc (n, file, name, index, read, scalar, offset)

      implicit none

      character (*) :: file, name
      character (120) :: fname, new_file_name, text

      integer i, index, iou, j, n, ib(10), ic(10)

      logical read, exists

      real offset, scalar, wnext, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "ctdbc.h"

      logical inqvardef

      real tmpij(imtm2,jmtm2)

      C2K = 273.15
      if (read) then
        fname = new_file_name (file)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Error => ", trim(fname), " does not exist."
          stop 'get_tdsbc in data.f'
        endif
        obc(:,:,n,inextm(n)) = c0
        ib(:) = 1
        ic(:) = 1
        ib(3) = inextd(n)
        ic(1) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        if (inqvardef(name, iou)) then
          call getvara (name, iou, imtm2*jmtm2, ib, ic, tmpij
     &,     scalar, offset)
          obc(2:imtm1,2:jmtm1,n,inextm(n)) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (obc(:,:,n,inextm(n)))
          text = "C"
          call getatttext (iou, name, 'units', text)
          if (name .ne. "L_diurtemp".and. trim(text) .eq. "K")
     &      obc(:,:,n,inextm(n)) = obc(:,:,n,inextm(n)) - C2K
          where (obc(:,:,n,inextm(n)).gt.1.e30) obc(:,:,n,inextm(n))=0.
        endif
      endif
      wnext = c1-wprev(n)
      do j=1,jmt
        do i=1,imt
          sbc(i,j,index) = wprev(n)*obc(i,j,n,iprevm(n))
     &                   + wnext*obc(i,j,n,inextm(n))
        enddo
      enddo
      n = n + 1

      return
      end
