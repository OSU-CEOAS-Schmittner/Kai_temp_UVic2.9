! source file: /usr/local/models/UVic_ESCM/2.9/updates/02/source/mom/checks.F
      subroutine checks (errorc, vmixset, hmixset)

      implicit none

      integer i, k, j, ip, kr, jq, n, num_mw, jb, jjs, jj, ncrow, jw
      integer je, in, is, numk

      real critv, t1, dymin, dxmin, jrow, dzmin, xlmax, dtxl, num
      real dxdymn, fimax, fmax, dysq, dxsq, clix, h1, h2, hx

      logical errorc, vmixset, hmixset
      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "hmixc.h"
      include "iounit.h"
      include "levind.h"
      include "isopyc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

!-----------------------------------------------------------------------
!     do consistency checks
!-----------------------------------------------------------------------

      write (stdout,'(/,20x,a,/)')
     &         'G E N E R A L    C O N S I S T E N C Y    C H E C K S'

      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "imt" less than 3 is not allowed         '
        errorc = .true.
      endif

      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "jmt" less than 4 is not allowed          '
        errorc = .true.
      endif

      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "consthmix"  cannot be enabled because another      '
     &,'            horizontal mixing scheme has been enabled          '
        errorc = .true.
      else
        hmixset = .true.
      endif
      if (.not.hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal mixing scheme has been enabled        '
        errorc = .true.
      endif

      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "constvmix"  cannot be enabled because another      '
     &,'            vertical mixing scheme has been enabled            '
        errorc = .true.
      else

!       set vmixset = true for enabeling "constvmix"

        vmixset = .true.
      endif

      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled          '
        errorc = .true.
      endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled          '
        errorc = .true.
      endif

      if (jmw .lt. 4) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 4 rows when using any'
     &,'           fourth order options                                '
        write (stdout,*)'           you have set jmw=',jmw
        errorc = .true.
      endif
      if (jmw .gt. 4) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 4 ("jmw"=4 will use the minimum memory)   '
        write (stdout,*)'             you have set jmw=',jmw
      endif
      if (jmw .gt. jmt) then
        write (stdout,9000)
     & '==> Error: the MW can not have more rows than "jmt"            '
        write (stdout,*)'           you have set jmw=',jmw, ', jmt=',jmt
        errorc = .true.
      endif
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Warning: The MW is open all the way ("jmw" = "jmt") which  '
     &,'             is the maximum memory configuration. Note that    '
     &,'             latitude rows are kept in the MW and not on disk! '
      endif

      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                         '
        write (stdout,*)'           nkflds is set = ',nkflds
        errorc = .true.
      endif

      if (dampts(1) .ne. c0 .or. dampts(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping time scale "dampts" is > zero but     '
     &,'             the "restorst" otpion is not enabled              '
      endif
      if (dampdz(1) .ne. c0 .or. dampdz(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping thickness "dampdz" is > zero but      '
     &,'             the "restorst" otpion is not enabled              '
      endif

      write (stdout,9000)
     & '==> Note: consthmix will only affect mixing of momentum        '
     &,'          since isopycmix was specified for tracer diffusion.  '
     &,'          kappa_h and Ah will be used as background mixing     '
     &,'          coefficients                                         '
      write (stdout,9000)
     & '==> Warning: although "sf_5_point" has no null space, it does  '
     &,'             not conserve total energy.                        '
      if ((ah+ahisop) .gt. 1.e11) then
        write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                  '
     &,'           "isopycmix" mixing option                           '
        errorc = .true.
      endif

      if (dtsf .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"      '
        errorc = .true.
      endif

      if (dtuv .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"      '
        errorc = .true.
      endif

      if (dtts .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"            '
        errorc = .true.
      endif

      critv = 1.e-6
      if (mod(rundays,dtts*secday) .gt. critv) then
        t1 = nint(rundays/(dtts*secday))*dtts*secday
        write (stdout,9000)
     & '==> Warning: there must be an integral number of density time  '
     &,'             steps within "rundays" (the integration time).    '
        write (stdout,*) '               (changed "rundays" from     '
     &,   rundays,' days to ', t1,' days to insure this condition)     '
          rundays = t1
      endif

      if (itmb) then
        write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in     '
     &,'             subsequent runs to prevent the time independent   '
     &,'             basin mask from being written more than once.     '
     &,'             This reduces the size of the diagnostic file.     '
      endif

      if (itrmb) then
        write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
      endif

      if (itavg) then
        write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
      endif
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "tmbint" is > zero but the '
     &,'             the "meridional_tracer_budget" option is not on.  '
      endif
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "xbtint"  is > zero but    '
     &,'             the "xbts" option is not enabled                  '
      endif
      if (dspint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "dspint"  is > zero but    '
     &,'             option "diagnostic_surf_height" is not enabled    '
      endif

      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient   '
     &,'             response is unimportant and multiple equilibria   '
     &,'             do not exist.                                     '
      endif

!     check for mixing coefficients larger than stability permits

      dymin  = dyt(2)
      dxmin  = dxt(2)
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
      enddo

      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the         '
     &,'             transient response is unimportant and multiple    '
     &,'             equilibria do not exist. stability tests will     '
     &,'             use "dtts" multiplied by the maximum "dtxcel"     '
      endif

      dtxl = dtts*xlmax
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
        fimax = 0.
        do k=1,km
          do i=2,imtm1
            fimax = max(fimax,fisop(i,j,k))
          enddo
        enddo
        if ((dtxl*(ah+ahisop*fimax))/dxdymn .ge. p25) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "ah" +    '
     &,'             "ahisop". use a smaller "dtts", "dtxcel", and/or  '
     &,'             "ah" + "ahisop"                                   '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
        if ((dtuv*am)/dxdymn .ge. p25) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "am".     '
     &,'             use a smaller "dtuv" and/or "am"                  '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
      if (dzt(1) .lt. 20.0e2) then
        write (stdout,9000)
     & '==> Warning: if shallow mixed layers develop, then enabling    '
     &,'             ifdef "shortwave" may help to deepen them. note   '
     &,'             that either you or the atmosphere must provide    '
     &,'             the solar short wave as a boundary condition.     '
      endif
      do k=1,km
        if ((dtts*dtxcel(k)*kappa_h)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_h" '
     &,'             use a smaller "dtts", "dtxcel", and/or "kappa_h"  '
       write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*kappa_m)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_m" '
     &,'             use a smaller "dtuv" and/or "kappa_m"             '
      endif
      write (stdout,9000)
     & '==> Warning: the full convective scheme is enabled.            '
     &,'             it will ignore "ncon" and remove all instability  '

!     check range of implicit factors

      if (acor .ne. 0) then
        write (stdout,9000)
     & '==> Error: "acor" must=0 when option damp_inertial_oscillation '
     &,'           is not enabled.                                     '
        errorc = .true.
      else

!       check for marginally resolved inertial oscillation

        fmax = epsln
        do jrow=2,jmtm1
          do i=2,imtm1
            fmax = max(fmax,abs(cori(i,jrow,1)))
          enddo
        enddo
        if (dtuv .gt. (1.0/6.0)*(c2*pi)/fmax) then
          write (stdout,9000)
     & '==> Error: the inertial oscillation is not resolved. reduce    '
     &,'           "dtuv" or use option "damp_inertial_oscillation"    '
          errorc = .true.
        endif
      endif

!-----------------------------------------------------------------------
!     search for topographic instabilities (based  on the  work of
!     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
!-----------------------------------------------------------------------

      num   = 50
      do j=2,jmtm1
        dysq = dyt(j)**2
        do i=2,imtm1
        if (kmu(i+1,j-1) .ne. 0 .and. kmu(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmu(i+1,j-1))
            h2   = zw(kmu(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write(stdout,*)
              write (stdout,'(a,a,i4,a,i4,a)')
     &        '==> Warning: Killworth topographic roughness condition'
     &,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     &         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
        enddo
      enddo

!     verify that the domain boundary is valid

      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,1) .ne. 0) is = i
        if (kmt(i,jmt) .ne. 0) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the southern boundary.                              '
        write (stdout,*) '           at j=1 and i=',is
      endif
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the northern boundary.                              '
        write (stdout,*) '           at j=jmt and i=',in
      endif

!     verify that each ocean point is at least 2 levels deep

      numk = 0
      do jrow=1,jmt
        do i=1,imt
          if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
            errorc = .true.
            write (stdout,*)
     &       ' Error: kmt(',i,',',jrow,') = 1 is not allowed    '
          endif
        enddo
      enddo
      if (numk .ne. 0) then
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean   '
     &,'           points.                                             '
      endif

      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
      if (errorc) stop '=>checks'

9000  format (/,(1x,a))

      return
      end
