! source file: /raid24/aschmitt/UVic2.9/MOBI1.9/nobio/updates/diago.F
      subroutine diago

!=======================================================================
!     write out diagnostic output
!=======================================================================

      implicit none

      integer numdsp, jrow, i, nislsp, mdscan, is, indp, ie, js, je
      integer io, m, mask, k, msk, num, n, iobadt, iobads, jte, jue
      integer jwte, jwue, ll, nv, ks, ke, iv, ih, lll, maxm, mloop
      integer ms, me, l, len, j, jj, iou

      logical defined

      real uext, vext, rnum, npt, dspcrt, rgrav, scl, reltim, rvolgk
      real rvolgt, rareag, avgper, fx, tnew, tsml, tbig, plicin, plicex
      real buoerr, enleak, dtconv, dtfilt, taux, tauy, erru, errv
      real contu, dchg, cont, tperiod, pwatts, csalt
      real tmbavg, convrt, sumsto, sumdiv, sumflx, sumdif, sumsor
      real sumvol, terror, serror, cmmday, cwatts, zmau, zmsmf,  zmsm
      real zmat, zmstf, zmst, tconv, tfilt, rnavg, timunit, tmp

      character(120) :: fname, ftbt, file_stamp, new_file_name
      character(32) :: nstamp
      save ftbt
      data ftbt  /' '/

      integer ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "ctmb.h"
      include "coord.h"
      include "cprnts.h"
      include "cregin.h"
      include "diag.h"
      include "docnam.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "stab.h"
      include "state.h"
      include "switch.h"
      include "tmngr.h"
      include "vmixc.h"
      include "cembm.h"
      include "ctavg.h"
      include "npzd.h"

!-----------------------------------------------------------------------
!     compute tracer averages under horizontal regions and write them
!     out
!-----------------------------------------------------------------------

      if (tavgts) then

!       initialize sums for horizontal regions & tracer averages

        do m=1,nt
          sumgf(m) = c0
          avggf(m) = c0
          sumgt(m) = c0
          avggt(m) = c0
          do mask=1,nhreg
            sumbt(mask,m) = c0
            avgbt(mask,m) = c0
            avgbf(mask,m) = c0
          enddo
          do k=1,km
            sumgk(k,m) = c0
            avggk(k,m) = c0
          enddo
        enddo

!       compute sums for tracer averages over horizontal regions

        do m=1,nt
          do mask=1,nhreg
            sumgf(m) = sumgf(m) + sumbf(mask,m)
            do k=1,km
              sumbt(mask,m) = sumbt(mask,m) + sumbk(mask,k,m)
              sumgk(k,m) = sumgk(k,m) + sumbk(mask,k,m)
            enddo
            sumgt(m) = sumgt(m) + sumbt(mask,m)
          enddo
        enddo

        do k=1,km
          if (volgk(k) .gt. c0) then
            rvolgk = c1 / volgk(k)
            do m=1,nt
              avggk(k,m) = sumgk(k,m) * rvolgk
              do mask=1,nhreg
                if (volbk(mask,k) .gt. c0) then
                  avgbk(mask,k,m) = sumbk(mask,k,m) / volbk(mask,k)
                endif
              enddo
            enddo
          endif
        enddo

        rvolgt = c1 / volgt
        rareag = c1 / areag
        do m=1,nt
          avggt(m) = sumgt(m) * rvolgt
          avggf(m) = sumgf(m) * rareag
          do mask=1,nhreg
            if (volbt(mask) .gt. c0) then
              avgbt(mask,m) = sumbt(mask,m) / volbt(mask)
            endif
            if (areab(mask) .gt. c0) then
              avgbf(mask,m) = sumbf(mask,m) / areab(mask)
            endif
          enddo
        enddo

!       write out regional tracer means

        if (iotavg .eq. stdout .or. iotavg .lt. 0) then
          write (stdout,'(//,30x,a,/)') 'T R A C E R    A V E R A G E S'
          do m=1,nt
            write(stdout,9004) trname(m), itt, stamp
            write(stdout,9001) (mask,mask=1,nhreg)
            do k=1,km
              write(stdout,9002) k, avggk(k,m),
     &                             (avgbk(mask,k,m),mask=1,nhreg)
            enddo
            write(stdout,9003) avggt(m), (avgbt(mask,m),mask=1,nhreg)
            write(stdout,9014) m, avggf(m), (avgbf(msk,m),msk=1,nhreg)
          enddo
        endif
        if (iotavg .ne. stdout .or. iotavg .lt. 0) then
          reltim = prelyr

          call getunit (io, 'tracer_avg.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Regional tracer averages written'
     &,      ' unformatted to file tracer_avg.dta on ts=',itt, ' ',stamp

          iotext = 'read (iotavg) reltim, nt, nhreg, km'
          write (io) pstamp, iotext, expnam
          write (io) reltim, nt, nhreg, km

          iotext = 'read (iotavg) ((avggk(k,n),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggk, km*nt)

          iotext =
     &    'read (iotavg) (((avgbk(l,k,n),l=1,nhreg),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbk, nhreg*km*nt)

          iotext = 'read (iotavg) (avggt(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggt, nt)

          iotext = 'read (iotavg) (avggf(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggf, nt)

          iotext = 'read (iotavg) ((avgbt(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbt, nhreg*nt)

          iotext = 'read (iotavg) ((avgbf(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbf, nhreg*nt)

          call relunit (io)
        endif
      endif
9001  format('    k','  All Regions ',9(1x,i7,5x))
9002  format(1x,i4,10(1x,e12.6))
9003  format('  AVG',10(1x,e12.6))
9004  format(/' Volume Weighted Averages for ',a12,' on ts =',i10,a33)
9014  format(/' FLX',i1,10(1x,e12.6),/)

!-----------------------------------------------------------------------
!     print integrals for monitoring the time step
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) call ta_mom_tsi (1)

      if (tsits .and. ntatio .gt. 0) then
        call ta_mom_tsi (2)

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (stdout,9601) itt, nstamp, tai_ek, tai_dt, tai_ds, tai_t
     &,                       tai_s, tai_tvar, tai_svar, tai_scan
        endif

        call def_tsi
        call def_tsi_mom (fname)
        tai_cocn = 0.
        tai_cfa2o = 0.

        avgper = tsiper*accel
        if (avgper .le. 1e-6) avgper = 0.
        tmp = 0.5
        time = time - tmp*avgper/365.

        call mom_tsi_out (fname, avgper, time, nstamp, tai_ek, tai_t
     &,                   tai_s, tai_tvar, tai_svar, tai_dt, tai_ds
     &,                   tai_scan, tai_hflx, tai_sflx, tai_dic
     &,                   tai_dicflx, tai_dic13, tai_dic13flx
     &,                   tai_alk, tai_o2, tai_o2flx
     &,                   tai_po4, tai_dop, tai_phyt, tai_zoop, tai_detr
     &,                   tai_no3, tai_don, tai_diaz, tai_din15
     &,                   tai_don15, tai_phytn15, tai_zoopn15
     &,                   tai_detrn15, tai_diazn15, tai_doc13
     &,                   tai_phytc13, tai_zoopc13, tai_detrc13
     &,                   tai_diazc13, tai_c14, tai_dc14
     &,                   tai_c14flx, tai_cfc11, tai_cfc11flx, tai_cfc12
     &,                   tai_cfc12flx, tai_otmax, tai_otmin
     &,                   tai_slh, tai_sspH, tai_ssCO3, tai_ssOc
     &,                   tai_ssOa, tai_sspCO2, tai_cocn, tai_cfa2o
     &,                   ntrec, tai_dfe, tai_ddfe
     &       )

        call ta_mom_tsi (0)
      endif
9601  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,' |dT|=',1pe13.6
     &,      ' |dS|=',1pe13.6,' Tbar=',1pe13.6,' Sbar=',1pe13.6
     &,      ' Tvar=',1pe13.6,' Svar=',1pe13.6, ' scans=',1pe13.6)

!-----------------------------------------------------------------------
!     show stability and CFL conditions, reynolds and peclet numbers
!-----------------------------------------------------------------------

      if (stabts) then
        write (stdout,'(///20x,a/15x,a,/)')
     & ' S T A B I L I T Y     A N A L Y S I S'
     &, '(The indicated locations are the most unstable ones)'

        write (stdout,*) ' longitudinal domain: ',cflons, ' to ',cflone
        write (stdout,*) ' latitudinal  domain: ',cflats, ' to ',cflate
        write (stdout,*) ' depth domain (m)   : ',cfldps*0.01
     &,                  ' to ',cfldpe*0.01

        write (stdout,'(/60x,a/)') ' CFL summary'
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local U velocity (',cflum,') is ',cflup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflu,',',jcflu,',',kcflu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflu),',',yu(jcflu),','
     &,  0.01*zt(kcflu),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local V velocity (',cflvm,') is ',cflvp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflv,',',jcflv,',',kcflv,')'
     &,  ' (lon,lat,dpt) = (',xu(icflv),',',yu(jcflv),','
     &,  0.01*zt(kcflv),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbu    (',cflwum,') is ',cflwup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwu,',',jcflwu,',',kcflwu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwu),',',yu(jcflwu),','
     &,  0.01*zw(kcflwu),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbt    (',cflwtm,') is ',cflwtp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwt,',',jcflwt,',',kcflwt,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwt),',',yu(jcflwt),','
     &,  0.01*zw(kcflwt),')'

        fx = 100.0
        if (cflup .gt. fx .or. cflvp .gt. fx .or. cflwup .gt. fx .or.
     &      cflwtp .gt. fx) then
          write (stdout,*)
     &      ' => Warning. CFL exceeded... computational mode exists!'
        endif

        write (stdout,'(/60x,a24/)') ' Reynolds number summary'

        write (stdout,10300) reynx, ireynx, jreynx, kreynx
     &,                      xu(ireynx), yu(jreynx), 0.01*zt(kreynx)
        write (stdout,10310) reynu, reynmu

        write (stdout,10400) reyny, ireyny, jreyny, kreyny
     &,                      xu(ireyny), yu(jreyny), 0.01*zt(kreyny)
        write (stdout,10410) reynv, reynmv

        write (stdout,10500) reynz, ireynz, jreynz, kreynz
     &,                      xu(ireynz), yu(jreynz), 0.01*zt(kreynz)
        write (stdout,10510) reynw, reynmw

        write (stdout,'(/60x,a22/)') ' Peclet number summary'

        write (stdout,10600) peclx, ipeclx, jpeclx, kpeclx
     &,                      xu(ipeclx), yu(jpeclx), 0.01*zt(kpeclx)
        write (stdout,10610) peclu, peclmu

        write (stdout,10700) pecly, ipecly, jpecly, kpecly
     &,                      xu(ipecly), yu(jpecly), 0.01*zt(kpecly)
        write (stdout,10710) peclv, peclmv

        write (stdout,10800) peclz, ipeclz, jpeclz, kpeclz
     &,                      xu(ipeclz), yu(jpeclz), 0.01*zt(kpeclz)
        write (stdout,10810) peclw, peclmw

!       show ficticious tracer extremums

        call getunit (iostab, 'iostab', 'formatted sequential rewind')
        rewind iostab
        write (stdout,'(/,10x,a/,11x,a,1pe10.3,a/)')
     &  'Spurious creation of local tracer extremum (if any) follow:'
     &, '(where tracer exceeds local extremum by ',tdig,'*tracer)'
        do num=1,1000
          read (iostab,'(i4, i4, i4, i2, 3g14.7)', end=101, err=101)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    't(i,k,jrow,n) = t(',i,',',k,',',jrow,',',n, ') = '
     &,   tnew,' : local min was ',tsml, ' : local max was ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 101
          endif
        enddo
101     continue
        call relunit (iostab)

!       show T and S outside allowable bounds used for density coeffs

        call getunit (iobadt, 'iobadt', 'formatted sequential rewind')
        rewind iobadt
        write (stdout,'(/,10x,a/)')
     &  'Temperatures (if any) outside allowable ranges follow:'
        do num=1,1000
          read (iobadt,'(i4, i4, i4, i2, 3g14.7)', end=201, err=201)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    't(i,k,jrow,n) = t(',i,',',k,',',jrow,',',n, ') = '
     &,   tnew,' : tmin is ',tsml, ' : tmax is ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 201
          endif
        enddo
201     continue
        call relunit (iobadt)

        call getunit (iobads, 'iobads', 'formatted sequential rewind')
        rewind iobads
        write (stdout,'(/,10x,a/)')
     &  'Salinities (if any) outside allowable ranges follow:'
        do num=1,1000
          read (iobads,'(i4, i4, i4, i2, 3g14.7)', end=301, err=301)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    'Error condition: t(i,k,jrow,n) = t(',i,',',k,',',jrow
     &,   ',',n, ') = ', tnew,' : smin is ',tsml, ' : smax is ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 301
          endif
        enddo
301     continue
        call relunit (iobads)
        write (stdout,'(/60x,a/)') ' End Stability Analysis'
      endif
10300 format (1x,'Maximim zonal Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10310 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10400 format (1x,'Maximim meridional Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10410 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10500 format (1x,'Maximim vertical Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10510 format (1x,' local Wu =',1pe9.2, ' and  mixing =',e9.2)
10600 format (1x,'Maximim zonal Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10610 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10700 format (1x,'Maximim meridional Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10710 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10800 format (1x,'Maximim vertical Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10810 format (1x,' local Wt =',1pe9.2, ' and  mixing =',e9.2)

!-----------------------------------------------------------------------
!     add external mode component of total work done
!-----------------------------------------------------------------------

      if (glents) then
        call ge3 (c2dtuv)
      endif

!-----------------------------------------------------------------------
!     integrate previously computed energy components vertically
!-----------------------------------------------------------------------

      if (glents) then
        jte  = 1
        jue  = 1
        jwte = 1
        jwue = 1
        do k=1,km
          wtlev(k,0) = c0
          wulev(k,0) = c0
        enddo
        do jrow=2,jmt-1
          do k=km,1,-1
            buoy(0,1) = buoy(0,1) + buoy(k,jrow)
            wtlev(k,0) = wtlev(k,0) + wtlev(k,jrow)
            wulev(k,0) = wulev(k,0) + wulev(k,jrow)
          enddo
          do ll=1,8
            do k=km,1,-1
              engint(0,ll,1) = engint(0,ll,1) + engint(k,ll,jrow)
            enddo
            engext(ll,1) = engext(ll,1) + engext(ll,jrow)
          enddo
          if (abs(tcerr(jrow)) .gt. abs(tcerr(jte))) jte = jrow
          if (abs(ucerr(jrow)) .gt. abs(ucerr(jte))) jue = jrow
          if (abs(wtbot(jrow)) .gt. abs(wtbot(jwte))) jwte = jrow
          if (abs(wubot(jrow)) .gt. abs(wubot(jwue))) jwue = jrow
        enddo
        buoy(0,1) = buoy(0,1)/ucellv
        do ll=1,8
          engint(0,ll,1) = engint(0,ll,1)/ucellv
          engext(ll,1)   = engext(ll,1)/ucellv
        enddo

        plicin = engint(0,1,1) - engint(0,2,1) - engint(0,3,1)
     &           - engint(0,4,1) - engint(0,5,1) - engint(0,6,1)
        plicex = engext(1,1) - engext(2,1) - engext(3,1)
     &            - engext(4,1)  - engext(5,1) - engext(6,1)
        buoerr = buoy(0,1) - engint(0,6,1) - engext(6,1)
        enleak = engint(0,2,1) + engint(0,3,1) + engext(2,1)
     &         + engext(3,1)

        if (ioglen .eq. stdout .or. ioglen .lt. 0) then
          write (stdout,'(///40x,a,/)') 'E N E R G Y    A N A L Y S I S'
          write (stdout,9100)
     &              'Globally averaged work done'
     &,                     ' for ts =', itt, stamp, ucellv, ucella(1)
          write (stdout,9101) ' time rate of change ',engint(0,1,1)
     &,                      engint(0,1,1), engext(1,1), engext(1,1)
          write (stdout,9101) ' horizontal advection', engint(0,2,1)
     &,                      engint(0,2,1), engext(2,1), engext(2,1)
          write (stdout,9101) ' vertical advection  ',engint(0,3,1)
     &,                      engint(0,3,1), engext(3,1), engext(3,1)
          write (stdout,9101) ' horizontal friction ',engint(0,4,1)
     &,                      engint(0,4,1), engext(4,1), engext(4,1)
          write (stdout,9101) ' vertical friction   ',engint(0,5,1)
     &,                      engint(0,5,1), engext(5,1), engext(5,1)
          write (stdout,9101) ' pressure forces     ',engint(0,6,1)
     &,                      engint(0,6,1), engext(6,1), engext(6,1)
          write (stdout,9101) ' ficticious sources  ',plicin
     &,                       plicin, plicex, plicex
          write (stdout,9101) ' work by wind        ',engint(0,7,1)
     &,                      engint(0,7,1), engext(7,1), engext(7,1)
          write (stdout,9101) ' bottom drag         ',engint(0,8,1)
     &,                      engint(0,8,1), engext(8,1), engext(8,1)
          write (stdout,9110) buoy(0,1), buoy(0,1), buoerr, buoerr
     &,                      enleak, enleak
          write (stdout,9111) tcerr(jte), itcerr(jte), jtcerr(jte)
     &,                       ktcerr(jte)
     &,                       ucerr(jue), iucerr(jue), jucerr(jue)
     &,                       kucerr(jue)
          write (stdout,9112) wtbot(jwte), iwtbot(jwte), jwtbot(jwte)
     &,                       kwtbot(jwte)
     &,                       wubot(jwue), iwubot(jwue), jwubot(jwue)
     &,                       kwubot(jwue)

          write (stdout,'(/a,a,//a,a,a,a,a,/1x,a,a)')
     &    'Average vertical velocity through bottom of "T" and "u" '
     &,   'cells at each level','level', '   adv_vbt err '
     &, '  "T" cell area ', '    adv_vbu err ','  "U" cell area'
     &, '(Note: adv_vbu err only goes to zero when non-zero values on'
     &, ' boundary cells are taken into account)'
          do k=1,km
            if (tcella(k) .ne. c0) then
              wtlev(k,0) = wtlev(k,0)/tcella(k)
            else
              wtlev(k,0) = c0
            endif
            if (ucella(k) .ne. c0) then
              wulev(k,0) = wulev(k,0)/ucella(k)
            else
              wulev(k,0) = c0
            endif
            write (stdout,'(i4,4(2x,e14.7))')
     &      k, wtlev(k,0), tcella(k), wulev(k,0), ucella(k)
          enddo
        endif

        if (ioglen .ne. stdout .or. ioglen .lt. 0) then
          reltim = prelyr

          call getunit (io, 'energy_int.dta'
     &,                'unformatted sequential append ieee')

          write (stdout, *)
     &      ' ==> Global energy integrals written unformatted'
     &,     ' to file energy_int.dta on ts =', itt, stamp
          iotext = 'read(ioglen) reltim,(engint(i),engext(i),i=1,8)'
          iotext(46:) = ',plicin,plicex,buoy,buoerr,enleak'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, (engint(0,i,1),engext(i,1),i=1,8)
     &,                 plicin, plicex, buoy(0,1), buoerr, enleak
          call relunit (io)
        endif
      endif
9100  format(///,1x,
     &/1x,a,a,i10,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', ocean surface area =',e16.9,' cm**2'//' work by:',14x
     &,     'internal mode                         external mode'/)
9101  format(a21,2(1pe15.6, ' (',z16,' hex)'))
9110  format(/' work by buoyancy forces   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' energy conversion error   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' nonlinear error           =',1pe14.6, ' (',z16,' hex)'/)
9111  format(/' max "t" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max "u" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')')
9112  format(/' max bottom "adv_vbt" (error)       =',1pe14.6
     &, ' at location', ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max bottom "adv_vbu" (slope vel)   =',1pe14.6
     &, ' at location',' (i,jrow,k) = ','(',i4,',',i4,',',i4,')')

!-----------------------------------------------------------------------
!     add the external mode part of d/dt into the momentum balance,
!     the external mode part of the implicit coriolis term, and the
!     surface pressure gradients into the specified volumes
!-----------------------------------------------------------------------

      if (trmbts) then

        call utb3

!-----------------------------------------------------------------------
!     integrate previously computed term balance components vertically
!-----------------------------------------------------------------------

        do n=0,numreg
          if (n .gt. 0) then
            nv = (n-1)/nhreg + 1
            ks = llvreg(nv,1)
            ke = llvreg(nv,2)
            do ll=1,17
              do k=ke,ks,-1
                termbm(0,ll,1,n) = termbm(0,ll,1,n) + termbm(k,ll,1,n)
                termbm(0,ll,2,n) = termbm(0,ll,2,n) + termbm(k,ll,2,n)
              enddo
            enddo
          else
            ks = 1
            ke = km
          endif
          do m=1,nt
            do ll=1,15
              do k=ke,ks,-1
                termbt(0,ll,m,n) = termbt(0,ll,m,n) + termbt(k,ll,m,n)
              enddo
            enddo

!           construct change due to convection and filtering

            dtconv = termbt(0,10,m,n) - termbt(0,9,m,n)
            dtfilt = termbt(0,1,m,n) - termbt(0,10,m,n)
            termbt(0,9,m,n) = dtconv
            termbt(0,10,m,n) = dtfilt
          enddo
        enddo

!       normalize integrals by appropriate volume (or area)

        do n=0,numreg
          do m=1,nt
            if (n .le. nhreg)  then
              stflx(m,n) = stflx(m,n)*rareat(n)
              asst(m,n)  = asst(m,n)*rareat(n)
            endif
            do ll=1,15
              termbt(0,ll,m,n) = termbt(0,ll,m,n)*rvolt(n)
            enddo
          enddo
        enddo

        do n=0,numreg
          if (n .le. nhreg) then
            smflx(1,0) = smflx(1,0) + smflx(1,n)
            smflx(2,0) = smflx(1,0) + smflx(2,n)
            smflx(1,n) = smflx(1,n)*rareau(n)
            smflx(2,n) = smflx(2,n)*rareau(n)
          endif
          if (n .gt. 0) then
            avgw(n) = avgw(n)*rvolu(n)
            do ll=1,17
              termbm(0,ll,1,n) = termbm(0,ll,1,n)*rvolu(n)
              termbm(0,ll,2,n) = termbm(0,ll,2,n)*rvolu(n)
            enddo
          endif
        enddo
        smflx(1,0) = smflx(1,0)*rareau(0)
        smflx(2,0) = smflx(2,0)*rareau(0)

        if (iotrmb .eq. stdout .or. iotrmb .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'T E R M    B A L A N C E S'

          n = 0
          taux = smflx(1,n)
          tauy = smflx(2,n)
          write (stdout,10106)
     &                   'All regions added together for ts ='
     &,                     itt, stamp, volt(n), areat(n)
     &,                     volu(n), areau(n)
          write (stdout,10104)
          write (stdout,10098) n, ' smf(1)   = ', taux,' dynes/cm**2   '
          write (stdout,10098) n, ' smf(2)   = ', tauy,' dynes/cm**2   '
          do m=1,nt
            write (stdout,10098) n, ustf(m,1), stflx(m,n), ustf(m,2)
          enddo
          write (stdout,10098) n, ' tot heat = '
     &,                           (termbt(0,15,1,n)*volt(n))
     &,                                        ' deg C * cm**3 '
          write (stdout,10098) n, ' sst      = ',asst(1,n)
     &,                                       ' deg C         '

          do n=0,numreg
            if (n .eq. 1) then
              write (stdout,10050)
     &  'Regional averaged Momentum & Tracer Term Balances for ts =   '
     &,       itt, stamp
            endif
            iv = 0
            if (n .gt. 0) then
              iv = (n-1)/nhreg + 1
              ih = n - (iv-1)*nhreg
              write (stdout,10100)
     &        'Momentum terms averaged over region #'
     &,          n, ': ', hregnm(ih), vregnm(iv), volu(n), areau(n)
              write (stdout,10104)
              write (stdout,10101) n,' dU/dt   = ', termbm(0,1,1,n)
     &,                              ' dV/dt   = ', termbm(0,1,2,n)
              write (stdout,10101) n,' -Px     = ', termbm(0,2,1,n)
     &,                              ' -Py     = ', termbm(0,2,2,n)
              write (stdout,10101) n,' -surf Px= ', termbm(0,12,1,n)
     &,                              ' -surf Py= ', termbm(0,12,2,n)
              write (stdout,10101) n,' -(UU)x  = ', termbm(0,3,1,n)
     &,                              ' -(UV)x  = ', termbm(0,3,2,n)
              write (stdout,10101) n,' -(VU)y  = ', termbm(0,4,1,n)
     &,                              ' -(VV)y  = ', termbm(0,4,2,n)
              write (stdout,10101) n,' -(WU)z  = ', termbm(0,5,1,n)
     &,                              ' -(WV)z  = ', termbm(0,5,2,n)
              write (stdout,10101) n,'ADV_Umet = ', termbm(0,13,1,n)
     &,                              '-ADV_Vmet= ', termbm(0,13,2,n)
              write (stdout,10101) n,'  DIFF_Ux= ', termbm(0,6,1,n)
     &,                              '  DIFF_Vx= ', termbm(0,6,2,n)
              write (stdout,10101) n,'  DIFF_Uy= ', termbm(0,7,1,n)
     &,                              '  DIFF_Vy= ', termbm(0,7,2,n)
              write (stdout,10101) n,'  DIFF_Uz= ', termbm(0,8,1,n)
     &,                              '  DIFF_Vz= ', termbm(0,8,2,n)
              write (stdout,10101) n,'DIFF_Umet= ', termbm(0,9,1,n)
     &,                              'DIFF_Vmet= ', termbm(0,9,2,n)
              write (stdout,10101) n,'  fV     = ', termbm(0,10,1,n)
     &,                              ' -fU     = ', termbm(0,10,2,n)
              write (stdout,10101) n,'  source = ', termbm(0,11,1,n)
     &,                              '  source = ', termbm(0,11,2,n)
              erru = c0
              errv = c0
              do lll=2,13
                erru = erru + termbm(0,lll,1,n)
                errv = errv + termbm(0,lll,2,n)
              enddo
              write (stdout,10101) n,'  error  = ', termbm(0,1,1,n)-erru
     &,                              '  error  = ', termbm(0,1,2,n)-errv

              write (stdout,*) ' '
              write (stdout,10101) n,' -U(U)x  = ', termbm(0,14,1,n)
     &,                              ' -U(V)x  = ', termbm(0,14,2,n)
              write (stdout,10101) n,' -V(U)y  = ', termbm(0,15,1,n)
     &,                              ' -V(V)y  = ', termbm(0,15,2,n)
              write (stdout,10101) n,' -W(U)z  = ', termbm(0,16,1,n)
     &,                              ' -W(V)z  = ', termbm(0,16,2,n)

!             mass conservation within volume: Ux + Vy + Wz

              contu = (termbm(0,3,1,n) + termbm(0,4,1,n) +
     &                termbm(0,5,1,n)) - (termbm(0,14,1,n) +
     &                termbm(0,15,1,n) + termbm(0,16,1,n))
              write (stdout,10101) n,' mass err= ', contu

              write (stdout,10101) n,'  ubar   = ', termbm(0,17,1,n)
     &,                              '  vbar   = ', termbm(0,17,2,n)
     &,                              '  wbar   = ', avgw(n)
            endif
            if (iv .eq. 1) then
              write (stdout,10101) n,'  surf Uz= ', smflx(1,n)
     &,                              '  surf Vz= ', smflx(2,n)
            endif

            if (n .eq. 0) then
              write (stdout,10051)
     &   'Global averaged (all basins) Tracer Term Balances for ts = '
     &,       itt, stamp
            else
              write (stdout,10100)
     &        'Tracer terms averaged over region   #'
     &,        n, ': ', hregnm(ih), vregnm(iv), volt(n), areat(n)
            endif

            do m=1,nt
              dchg            = termbt(0,2,m,n) + termbt(0,3,m,n) +
     &                          termbt(0,4,m,n) + termbt(0,5,m,n) +
     &                          termbt(0,6,m,n) + termbt(0,7,m,n) +
     &                          termbt(0,8,m,n) + termbt(0,9,m,n) +
     &                          termbt(0,10,m,n)
              terr(m)        = termbt(0,1,m,n) - dchg
            enddo
            maxm = (nt-1)/7 + 1
            do mloop=1,maxm
              ms = (mloop-1)*7 + 1
              me = min(ms + 6,nt)
              write (stdout,10103) (trname(m),m=ms,me)
              write (stdout,10102) n,' dT/dt   = ', (termbt(0,1,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(UT)x  = ', (termbt(0,2,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(VT)y  = ', (termbt(0,3,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(WT)z  = ', (termbt(0,4,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tx= ', (termbt(0,5,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Ty= ', (termbt(0,6,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tz= ', (termbt(0,7,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  source = ', (termbt(0,8,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  convct = ', (termbt(0,9,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  filter = ', (termbt(0,10,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  error  = ', (terr(m)
     &,                                           m=ms,me)

              write (stdout,*) ' '
              write (stdout,10102) n,' -U(T)x  = ', (termbt(0,11,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -V(T)y  = ', (termbt(0,12,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -W(T)z  = ', (termbt(0,13,m,n)
     &,                                           m=ms,me)

!             mass conservation within volume: Ux + Vy + Wz

              cont = (termbt(0,2,1,n) + termbt(0,3,1,n) +
     &               termbt(0,4,1,n)) - (termbt(0,11,1,n) +
     &               termbt(0,12,1,n) + termbt(0,13,1,n))
              write (stdout,10102) n,' mass err= ', cont

              write (stdout,10102) n,'  chg var= ', (termbt(0,14,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  tbar   = ', (termbt(0,15,m,n)
     &,                                           m=ms,me)
              if (iv .eq. 1) then
                write (stdout,10102) n,'  surflx = ', (stflx(m,n)
     &,                                           m=ms,me)
                taux = smflx(1,n)
                tauy = smflx(2,n)
                write (stdout,10105) ' Regionally averaged quantities:'
                write (stdout,'(1x)')
                write (stdout,10098) n
     &,            ' smf(1)   = ', taux,' dynes/cm**2   '
                write (stdout,10098) n
     &,            ' smf(2)   = ', tauy,' dynes/cm**2   '
                do m=1,nt
                  write (stdout,10098) n
     &,            ustf(m,1), stflx(m,n), ustf(m,2)
                enddo
                if (ms .eq. 1) then
                  write (stdout,10098) n, ' tot heat = '
     &,             (termbt(0,15,1,n)*volt(n)),' deg C * cm**3 '
                  write (stdout,10098) n
     &,            ' sst      = ', asst(1,n),' deg C         '
                endif
              endif
            enddo
          enddo
        endif

        if (iotrmb .ne. stdout .or. iotrmb .lt. 0) then

          reltim  = prelyr
          tperiod = 0.0
          call getunit (io, 'term_bal.dta'
     &,                'unformatted sequential append ieee')

          write (stdout, *)
     &     ' ==> Term balances written unformatted to file term_bal.dta'
     &,    ' for ts =', itt, ' ',stamp

          iotext =' read (iotrmb) reltim, nt, numreg, nhreg, km'
          iotext(45:) =', ((ustf*15(n,i),n=1,nt),i=1,2)'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, nt, numreg, nhreg, km, ustf

          iotext ='read (iotrmb) (hregnm*40(n),n=1,nhreg)'
          iotext(39:)=', (vregnm*20(n),n=1,nvreg)'
          write (io) pstamp, iotext, expnam
          write (io)  hregnm, vregnm

          iotext ='read (iotrmb) (trname*12(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          write (io)  trname

          iotext ='read (iotrmb) (volu(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volu, numreg+1)

          iotext ='read (iotrmb) (volt(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volt, numreg+1)

          iotext ='read (iotrmb) (areau(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areau, numreg+1)

          iotext ='read (iotrmb) (areat(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areat, numreg+1)

          iotext ='read (iotrmb) ((smflx(i,l),i=1,2),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, smflx, 2*(nhreg+1))

          iotext ='read (iotrmb) ((stflx(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, stflx, nt*(nhreg+1))

          iotext ='read (iotrmb) ((asst(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, asst, nt*(nhreg+1))

          iotext ='read (iotrmb) (avgw(l),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgw, numreg)

          iotext = 'read (iotrmb) (((termbm'
          iotext(24:) = '(i,n,l),i=1,17),n=1,2),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbm(0,i,n,l),i=1,17),n=1,2),l=1,numreg)

          iotext = 'read (iotrmb) (((termbt'
          iotext(24:) = '(i,n,l),i=1,15),n=1,nt),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbt(0,i,n,l),i=1,15),n=1,nt),l=0,numreg)

          call relunit (io)
        endif
      endif
10050 format(///,1x,/1x,a,i10,a/)
10051 format(///,1x,a,i10,a/)
10098 format(1x,'(',i4,')',a12,(1pe16.8,a15))
10100 format(/1x,a,i5,a,a,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', horizontal ocean area =',e16.9,' cm**2'/)
10101 format(1x,'(',i4,')',a11,1pe15.7, 2x, a11,1pe15.7, 2x, a11
     &,      1pe15.7)
10102 format(1x,'(',i4,')',a11,7(1pe15.7))
10103 format(' Region',11x,8a15)
10104 format(' Region')
10105 format(/8x,a32)
10106 format(///,1x,
     &/1x,a,i10,a/1x,'ocean "t" volume =',e16.9,' cm**3 '
     &, ', ocean "t" surface area =',e16.9,' cm**2 '
     &/ 1x,'ocean "u" volume =',e16.9,' cm**3, '
     &, ' ocean "u" surface area =',e16.9,' cm**2'/)

!-----------------------------------------------------------------------
!     write out the gyre_components
!-----------------------------------------------------------------------

!     convert heat transport to petawatts,
!     salt transport to 10**10 cm**3/sec

      if (gyrets) then
        pwatts = 4.186e-15
        csalt  = 1.e-10
        do jrow=1,jmt
          do ll=1,8
            ttn(ll,jrow,1)=ttn(ll,jrow,1)*pwatts
            ttn(ll,jrow,2)=ttn(ll,jrow,2)*csalt
          enddo
        enddo
        if (iogyre .eq. stdout .or. iogyre .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'G Y R E   C O M P O N E N T S'

          write (stdout,8195)
          do jrow=2,jmtm2
            l = jmt - jrow
            write (stdout,8196) l, (ttn(i,l,1),i=1,8)
     &,                            (ttn(i,l,2),i=1,8)
          enddo

          do m=0,nhreg
            do jrow=1,jmt
              do ll=6,8
                ttn2(ll,jrow,1,m) = ttn2(ll,jrow,1,m)*4.186e-15
                ttn2(ll,jrow,2,m) = ttn2(ll,jrow,2,m)*1.e-10
              enddo
            enddo
          enddo

          write (stdout,82001)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
              write (stdout,8204) jrow,(ttn2(i,jrow,1,m),i=6,8)
            enddo
          enddo

          write (stdout,8205)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
              write (stdout,8204) jrow,(ttn2(i,jrow,2,m),i=6,8)
            enddo
          enddo
        endif
        if (iogyre .ne. stdout .or. iogyre .lt. 0) then
          reltim = prelyr

          call getunit (io, 'gyre_comp.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Gyre components written'
     &,      ' unformatted to file gyre_comp.dta on ts=',itt, ' ',stamp

          iotext = 'read (iogyre) jmt, ntmin2, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, ntmin2, reltim

          iotext =
     &      'read (iogyre) (((ttn(l,j,n),l=1,8),j=1,jmt),n=1,ntmin2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, ttn, 8*jmt*ntmin2)

          iotext =
     &    'read (iogyre) (((ttn2(6..8,1..jmt,1..2,0..nhreg)'
          write (io) pstamp, iotext, expnam
          len = 3*jmt*nt*(1+nhreg)
          call wrufio (io, ttn2, len)

          call relunit (io)
        endif
      endif
8195  format(//,30x,'Gyre Components'/
     &,6x,'northward transport of heat (x10**15 watts)'
     &       ,22x,'northward transport of salt (x10**10 cm**3/sec)',/,
     &       3x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     &       'diffus   total'))
8196  format(i4,8f8.3,1x,8f8.3)
82001 format (/,6x,'northward transport of heat (x10**15 watts)',/)
8201  format (/,22x,a40,/)
8202  format (/,22x,' Global ',/)
8203  format (8x,'total advection',2x,'total diffusion',2x,'total')
8204  format (2x,i3,3x,3(e12.6,5x))
8205  format (/,6x,'northward transport of salt (x10**10 cm**3/sec)',/)

!-----------------------------------------------------------------------
!     write out the meridional overturning (mass)
!-----------------------------------------------------------------------

      if (vmsfts) then
        if (iovmsf .eq. stdout .or. iovmsf .lt. 0) then

          write (stdout,'(///,40x,a,/)')
     &    'M E R I D I O N A L    O V E R T U R N I N G'

          scl = 1.e12
          write (stdout,8194)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
          call matrix (vmsf, jmt, js, je, 1, km, scl)
        endif
        if (iovmsf .ne. stdout .or. iovmsf .lt. 0) then
          reltim = prelyr

          call getunit (io, 'overturn.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Meridional overturning of mass'
     &,    ' written unformatted to file overturn.dta on ts=',itt
     &,    ' ',stamp

          iotext = 'read (iovmsf) jmt, km, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, km, reltim

          iotext = 'read (iovmsf)  (zw(k),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zw, km)

          iotext = 'read (iovmsf)  (yu(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yu, jmt)

          iotext = 'read (iovmsf)  ((vmsf(j,k),j=1,jmt),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, vmsf, jmt*km)

          call relunit (io)
        endif
      endif
8194  format(/,' meridional overturning stream function (sverdrups)')

!-----------------------------------------------------------------------
!     save the time mean averages on the "averaging" grid
!-----------------------------------------------------------------------

      if (timavgts .and. timavgint .gt. c0) then
        call avgout
      endif

      return
      end

      subroutine ta_mom_tsi (m)

!=======================================================================
!     ocean data time integral averaging

!     input:
!       m = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer i, j, k, ib(10), ic(10), m, n

      real dens, tq, sq, drodt, drods, drhodt, drhods, ddensdtdt
      real ddensdtds, ddensdsds, rntatio, rnv_otsf, rnt_slh, tarea
      real cstdyt, area, s, tins, tinsit, dins, sum, d, tmp, tmp1, tmp2

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "coord.h"
      include "cregin.h"
      include "grdvar.h"
      include "diag.h"
      include "diaga.h"
      include "emode.h"
      include "levind.h"
      include "mw.h"
      include "iounit.h"
      include "state.h"
      include "dens.h"

      real otsf(km)
      real dmsk(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

      if (m .eq. 0.) then

!       zero
        ntatio = 0
        tai_otmax = 0.
        tai_otmin = 0.
        tai_slh = 0.
        tai_ek = 0.
        tai_t = 0.
        tai_s = 0.
        tai_tvar = 0.
        tai_svar = 0.
        tai_dt = 0.
        tai_ds = 0.
        tai_scan = 0.
        tai_hflx = 0.
        tai_sflx = 0.
        tai_dic = 0.
        tai_dicflx = 0.
        tai_dic13 = 0.
        tai_dic13flx = 0.
        tai_dicwflx = 0.
        tai_dic13wflx = 0.
        tai_alk = 0.
        tai_o2 = 0.
        tai_o2flx = 0.
        tai_po4 = 0.
        tai_dop = 0.
        tai_phyt = 0.
        tai_zoop = 0.
        tai_detr = 0.
        tai_no3 = 0.
        tai_dfe = 0.
        tai_ddfe = 0.
        tai_don = 0.
        tai_diaz = 0.
        tai_din15 = 0.
        tai_don15 = 0.
        tai_phytn15 = 0.
        tai_coccn15 = 0.
        tai_zoopn15 = 0.
        tai_detrn15 = 0.
        tai_diazn15 = 0.
        tai_doc13 = 0.
        tai_phytc13 = 0.
        tai_coccc13 = 0.
        tai_zoopc13 = 0.
        tai_detrc13 = 0.
        tai_diazc13 = 0.
        tai_c14 = 0.
        tai_dc14 = 0.
        tai_c14flx = 0.
        tai_sspH = 0.
        tai_ssCO3 = 0.
        tai_ssOc = 0.
        tai_ssOa = 0.
        tai_sspCO2 = 0.
        tai_c = 0.
        tai_caco3 = 0.
        tai_caco3c13 = 0.
        tai_d_B = 0.
      elseif (m .eq. 1) then

!       accumulate
        ntatio = ntatio + 1

        tmp1 = 0.
        tmp2 = 0.
        if (nv_otsf .gt. 0) then
          tmp1 = -1.e20
          tmp2 = 1.e20
          rnv_otsf = 1./float(nv_otsf)
          do j=jsot,jeot
            do k=1,km
              v_otsf(j,k) = v_otsf(j,k)*rnv_otsf*dzt(k)*csu(j)
            enddo
!           integrate up from the bottom
            otsf(km) = 0.
            do k=km-1,1,-1
              otsf(k) = otsf(k+1) - v_otsf(j,k)
            enddo
            do k=ksot,keot
              if (otsf(k) .gt. tmp1) tmp1 = otsf(k)
              if (otsf(k) .lt. tmp2) tmp2 = otsf(k)
            enddo
          enddo
        endif
        tai_otmax = tai_otmax + tmp1
        tai_otmin = tai_otmin + tmp2
        nv_otsf = 0
        v_otsf(:,:) = 0.

        tmp = 0.
        if (nt_slh .gt. 0) then
          t_slh(:,:,:,:) = t_slh(:,:,:,:)/float(nt_slh)
          tarea = 0.
          do j=2,jmtm1
            cstdyt = cst(j)*dyt(j)
            do i=2,imtm1
              if (kmt(i,j) .ge. 1) then
                area = dxt(i)*cstdyt
                tarea = tarea + area
                do k=1,kmt(i,j)
                  s = (t_slh(i,j,k,2)*1000.) + 35.
                  d = zt(k)*0.01
                  tins = tinsit(t_slh(i,j,k,1), s, d)
                  dins = dens(tins-to(k), t_slh(i,j,k,2)-so(k), k)
!                 volume integrate (rho - rho0)/rho0
!                 assume rho0 = 1g/cm3 so difference is just delta rho
!                 dslh = sum[((rho1-rho0)/rho0)*dvol]/area
!                      - sum[((rho2-rho0)/rho0)*dvol}/area
!                      = sum[((rho1-rho2)/rho0)*dvol]/area
                  tmp = tmp - (dins - d_slh(i,j,k))*area*dzt(k)
                enddo
              endif
            enddo
          enddo
          if (tarea .gt. 0) tmp = tmp/tarea
        endif
        tai_slh = tai_slh + tmp
        nt_slh = 0
        t_slh(:,:,:,:) = 0.

        do j=2,jmtm1
          do k=km,1,-1
            ektot(0,1) = ektot(0,1) + ektot(k,j)
          enddo
          do n=1,nt
            do k=km,1,-1
              tbar(0,n,1)   = tbar(0,n,1) + tbar(k,n,j)
              travar(0,n,1) = travar(0,n,1) + travar(k,n,j)
              dtabs(0,n,1)  = dtabs(0,n,1) + dtabs(k,n,j)
            enddo
          enddo
        enddo
        ektot(0,1) = ektot(0,1)/ucellv

        do n=1,nt
          tbar(0,n,1) = tbar(0,n,1)/tcellv
          travar(0,n,1) = travar(0,n,1)/tcellv - tbar(0,n,1)**2
          dtabs(0,n,1) = dtabs(0,n,1)/tcellv
        enddo

        tai_ek = tai_ek + ektot(0,1)
        tai_t = tai_t + tbar(0,itemp,1)
        tai_s = tai_s + tbar(0,isalt,1)
        tai_tvar = tai_tvar + travar(0,itemp,1)
        tai_svar = tai_svar + travar(0,isalt,1)
        tai_dt = tai_dt + dtabs(0,itemp,1)
        tai_ds = tai_ds + dtabs(0,isalt,1)
        tai_scan = tai_scan + mscan

        dmsk(:,:) = 1.
        where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.

        call areaavg (sbc(1,1,ihflx), dmsk, tmp)
        tai_hflx = tai_hflx + tmp
        call areaavg (sbc(1,1,isflx), dmsk, tmp)
        tai_sflx = tai_sflx + tmp

      elseif (m .eq. 2 .and. ntatio .ne. 0) then

!       average
        rntatio = 1./float(ntatio)
        tai_otmax = tai_otmax*rntatio
        tai_otmin = tai_otmin*rntatio
        tai_slh = tai_slh*rntatio
        tai_ek = tai_ek*rntatio
        tai_t = tai_t*rntatio
        tai_s = tai_s*rntatio
        tai_tvar = tai_tvar*rntatio
        tai_svar = tai_svar*rntatio
        tai_dt = tai_dt*rntatio
        tai_ds = tai_ds*rntatio
        tai_scan = tai_scan*rntatio
        tai_hflx = tai_hflx*rntatio
        tai_sflx = tai_sflx*rntatio
        tai_dic = tai_dic*rntatio
        tai_dicflx = tai_dicflx*rntatio
        tai_dic13 = tai_dic13*rntatio
        tai_dic13flx = tai_dic13flx*rntatio
        tai_dicwflx = tai_dicwflx*rntatio
        tai_dic13wflx = tai_dic13wflx*rntatio
        tai_alk = tai_alk*rntatio
        tai_o2 = tai_o2*rntatio
        tai_o2flx = tai_o2flx*rntatio
        tai_po4 = tai_po4*rntatio
        tai_dop = tai_dop*rntatio
        tai_phyt = tai_phyt*rntatio
        tai_zoop = tai_zoop*rntatio
        tai_detr = tai_detr*rntatio
        tai_no3 = tai_no3*rntatio
        tai_don = tai_don*rntatio
        tai_diaz = tai_diaz*rntatio
        tai_din15 = tai_din15*rntatio
        tai_don15 = tai_don15*rntatio
        tai_phytn15 = tai_phytn15*rntatio
        tai_coccn15 = tai_coccn15*rntatio
        tai_zoopn15 = tai_zoopn15*rntatio
        tai_detrn15 = tai_detrn15*rntatio
        tai_diazn15 = tai_diazn15*rntatio
        tai_doc13 = tai_doc13*rntatio
        tai_phytc13 = tai_phytc13*rntatio
        tai_zoopc13 = tai_zoopc13*rntatio
        tai_detrc13 = tai_detrc13*rntatio
        tai_diazc13 = tai_diazc13*rntatio
        tai_c14 = tai_c14*rntatio
        tai_dc14 = tai_dc14*rntatio
        tai_c14flx = tai_c14flx*rntatio
        tai_cfc11 = tai_cfc11*rntatio
        tai_cfc11flx = tai_cfc11flx*rntatio
        tai_cfc12 = tai_cfc12*rntatio
        tai_cfc12flx = tai_cfc12flx*rntatio
        tai_sspH = tai_sspH*rntatio
        tai_ssCO3 = tai_ssCO3*rntatio
        tai_ssOc = tai_ssOc*rntatio
        tai_ssOa = tai_ssOa*rntatio
        tai_sspCO2 = tai_sspCO2*rntatio
        tai_dfe = tai_dfe*rntatio
        tai_ddfe = tai_ddfe*rntatio
      endif

      return
      end
      real function tinsit (ti, si, di)
!=======================================================================
!     function to calculate insitu temperature using Newtons method

!     input:
!       ti = potential temperature (C)
!       si = salinity (psu)
!       di = depth (m)

!     output:
!       tinsit = institu temperataure
!=======================================================================

      implicit none

      integer iter

      real, intent(in) :: di, si, ti
      real(kind=8) d, dfdt, error, s, t, ta, tb, tol

      error = 1.
      iter = 0
      tol = 1.e-5
      t = ti
      d = di
      s= si
      do while (abs(error) .gt. tol .and. iter .lt. 25)
        iter = iter + 1
        call potem(t, s, d, ta)
        call potem(t+1.e-3 ,s, d, tb)
        dfdt = (tb - ta)*1.e3
        error = (ta - ti)/dfdt
        t = t - error
      enddo
      if (abs(error) .gt. tol) then
        print*, "=> Error: insitu temperature did not converge"
        print*, "          ti = ", ti
        print*, "          si = ", si
        print*, "          di = ", di
        print*, "          error = ", abs(error)
        print*, "          tinsit = ", t
        print*, "   Stopped in function tinsit in diago.F"
        stop
      endif
      tinsit = t

      return
      end
