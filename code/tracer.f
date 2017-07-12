! source file: /raid24/aschmitt/UVic2.9/karin/mwc15_npzd_fe_n15_c13_alk_caco3/updates/tracer.F

      subroutine tracer (joff, js, je, is, ie)

!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer istrt, iend, i, k, j, ip, kr, jq, n, jp, jrow, iou, js
      integer je, limit, joff, is, ie, kmx, m, kb, idiag, index
      integer it(10), iu(10), ib(10), ic(10), nfnpzd, mfnpzd, mxfnpzd
      integer id_time, id_xt, id_yt, id_zt,fe_jlo,fe_m,fe_k,fe_n

      logical inqvardef, exists

      parameter (fe_n = 14)

      real time
      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, zmax, cont, drho
      real drhom1, wt, ahbi_cstr, ahbi_csu_dyur, gamma, fy, fyz
      real fe_dy, fe_conc, fe_x(fe_n), fe_y(fe_n), bctz

      real expofe, impofe, feorgads, remife, thetamax, deffe, fecol
      real thetachl, chl, chl_D, feprime, fesed, bfe, sgb, thetamax_C
      real deffe_C

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "timeavgs.h"
      include "tmngr.h"
      include "vmixc.h"
      include "diaga.h"
      real rrc14std
      include "ice.h"
      include "npzd.h"
      real t_in(km), po4_in(km)
      real o2_in(km)
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
      real sgb_in(km)
      real no3_in(km)
      real din15_in(km)
      real expo, tnpzd(km,ntnpzd)
      include "atm.h"
      include "cembm.h"

      parameter (istrt=2, iend=imt-1)

      real twodt(km)
      real ai, hi, hs, rctheta, dayfrac, swr, declin
      real src(imt,km,jsmw:jemw,nsrc)
      include "isopyc.h"
      include "fdift.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
        enddo
      enddo

!-----------------------------------------------------------------------
!     ocean biogeochemistry and isotopes
!-----------------------------------------------------------------------
      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      nbio = c2dtts/dtnpzd
      dtbio = c2dtts/nbio
      rdtts = 1./c2dtts
      rnbio = 1./nbio

      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            ai = aice(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsno(i,jrow,2)
!     calculate day fraction and incoming solar
!     angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
!     KK check this
            swr = (tap+tap_C)*dnswr(i,jrow)*1e-3*(1. + ai*(exp(-ki
     &          *(hi + hs)) - 1.))

!        set prognostic variables
            tnpzd(:,imobipo4) = t(i,:,j,ipo4,taum1)
            tnpzd(:,imobiphyt) = t(i,:,j,iphyt,taum1)
            tnpzd(:,imobizoop) = t(i,:,j,izoop,taum1)
            tnpzd(:,imobidetr) = t(i,:,j,idetr,taum1)
            tnpzd(:,imobidic) = t(i,:,j,idic,taum1)
            tnpzd(:,imobidop) = t(i,:,j,idop,taum1)
            tnpzd(:,imobino3) = t(i,:,j,ino3,taum1)
            tnpzd(:,imobidon) = t(i,:,j,idon,taum1)
            tnpzd(:,imobidiaz) = t(i,:,j,idiaz,taum1)
            tnpzd(:,imobidin15) = t(i,:,j,idin15,taum1)
            tnpzd(:,imobidon15) = t(i,:,j,idon15,taum1)
            tnpzd(:,imobiphytn15) = t(i,:,j,iphytn15,taum1)
            tnpzd(:,imobizoopn15) = t(i,:,j,izoopn15,taum1)
            tnpzd(:,imobidetrn15) = t(i,:,j,idetrn15,taum1)
            tnpzd(:,imobidiazn15) = t(i,:,j,idiazn15,taum1)
            tnpzd(:,imobicoccn15) = t(i,:,j,icoccn15,taum1)
            tnpzd(:,imobicocc) = t(i,:,j,icocc,taum1)
            tnpzd(:,imobicaco3) = t(i,:,j,icaco3,taum1)
            tnpzd(:,imobidfe) = t(i,:,j,idfe,taum1)
            tnpzd(:,imobidetrfe) = t(i,:,j,idetrfe,taum1)
            tnpzd(:,imobidic13) = t(i,:,j,idic13,taum1)
            tnpzd(:,imobidoc13) = t(i,:,j,idoc13,taum1)
            tnpzd(:,imobiphytc13) = t(i,:,j,iphytc13,taum1)
            tnpzd(:,imobizoopc13) = t(i,:,j,izoopc13,taum1)
            tnpzd(:,imobidetrc13) = t(i,:,j,idetrc13,taum1)
            tnpzd(:,imobidiazc13) = t(i,:,j,idiazc13,taum1)
            tnpzd(:,imobicoccc13) = t(i,:,j,icoccc13,taum1)
            tnpzd(:,imobicaco3c13) = t(i,:,j,icaco3c13,taum1)
!        set other input variables
            t_in(:) = t(i,:,j,itemp,taum1) !degree Celsius
            po4_in(:) = t(i,:,j,ipo4,taum1)
            o2_in(:) = t(i,:,j,io2,taum1)*1000. ! uM
            s_in(:) = 1.e3*t(i,:,j,isalt,taum1) + 35.0
            dic_in(:) = t(i,:,j,idic,taum1) !mol/m^3
            alk_in(:)  = t(i,:,j,ialk,taum1) !eq/m^3
            co2_in = co2ccn
            dic13_in(:) = tnpzd(:,imobidic13) !mol/m^3
            no3_in(:) = tnpzd(:,imobino3)
            din15_in(:) = tnpzd(:,imobidin15)
            sgb_in(:) = sg_bathy(i,j,:)
!  end input
            call mobi_driver(
     &                 kmt(i,jrow), c2dtts, rctheta, dayfrac, swr, tnpzd
     &,                t_in, po4_in
     &,                o2_in
     &,                s_in, dic_in, alk_in, co2_in
     &,                dic13_in
     &,                sgb_in
     &,                no3_in
     &,                din15_in
     &,                src(i,:,j,:))

!juan added dust iron here instead of in the set source/sinks part:
! comes in mol/m2/s gets converted to mmol/m3/s
                src(i,1,j,isdfe) = src(i,1,j,isdfe) +
     &          sbc(i,j,idfeadep)*1000./(dzt(1)/100.)  !dzt is in cm
!Add hydrothermal Fe:
               do k=1,kmt(i,j)
                src(i,k,j,isdfe) = src(i,k,j,isdfe) +
     &          fe_hydr(i,j,k)
               enddo

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------
            if (timavgperts .and. .not. euler2) then
             expo = rprca
             do k = 1, kmt(i,jrow)
              ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(k)
              ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(k)
              ta_rgraz_Z(i,k,jrow) = ta_rgraz_Z(i,k,jrow) + rgraz_Z(k)
              ta_rgraz_Det(i,k,jrow) = ta_rgraz_Det(i,k,jrow)
     &                                 + rgraz_Det(k)
              ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(k)
              ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(k)
              ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(k)
              ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(k)
              ta_rnpp_dop(i,k,jrow) = ta_rnpp_dop(i,k,jrow)
     &             + rnpp_dop(k)
              ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow) + rnpp_D(k)
              ta_rnpp_D_dop(i,k,jrow) = ta_rnpp_D_dop(i,k,jrow)
     &             + rnpp_D_dop(k)
              ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow) + rgraz_D(k)
              ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow) + rmorp_D(k)
              ta_rmorpt_D(i,k,jrow) = ta_rmorpt_D(i,k,jrow)
     &                               + rmorpt_D(k)
              ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(k)
              ta_rbdeni(i,k,jrow) = ta_rbdeni(i,k,jrow) + rbdeni(k)
              ta_rnpp_C_dop(i,k,jrow) = ta_rnpp_C_dop(i,k,jrow)
     &             + rnpp_C_dop(k)
!       print*,'KK tracer', rcalatt(k), rcalpro(k)
                ta_rcalatt(i,k,jrow) = ta_rcalatt(i,k,jrow)
     &                               + rcalatt(k)
                ta_rprocal(i,k,jrow) = ta_rprocal(i,k,jrow)
     &                               + rcalpro(k)
                ta_rnpp_C(i,k,jrow) = ta_rnpp_C(i,k,jrow)
     &                              + rnpp_C(k)
                ta_rgraz_C(i,k,jrow) = ta_rgraz_C(i,k,jrow)
     &                               + rgraz_C(k)
                ta_rmorp_C(i,k,jrow) = ta_rmorp_C(i,k,jrow)
     &                               + rmorp_C(k)
                ta_rmorpt_C(i,k,jrow)= ta_rmorpt_C(i,k,jrow)
     &                               + rmorpt_C(k)
              ta_rdissl(i,k,jrow) = ta_rdissl(i,k,jrow)
     &                              + rdissl(k)
              ta_rexpocaco3(i,k,jrow) =
     &          ta_rexpocaco3(i,k,jrow) + rexpocaco3(k)
              expo = expo*dztr(k)
              ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(k)
              ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(k)
              expo = expo - rprca*rcak(k)
              ta_rexpocal(i,k,jrow) = ta_rexpocal(i,k,jrow) + expo
              expo = expo*dzt(k)
              ta_rwcdeni(i,k,jrow) = ta_rwcdeni(i,k,jrow)
     &           + rwcdeni(k)
                ta_rremife(i,k,jrow) = ta_rremife(i,k,jrow)
     &               + rremife(k)
                ta_rexpofe(i,k,jrow) = ta_rexpofe(i,k,jrow)
     &               + rexpofe(k)
!
             enddo
            endif
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt

!-----------------------------------------------------------------------
!       calculate advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
     &                         ah_cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
     &                         diff_cnt*
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
            diff_fb(i,0,j)  = stf(i,j,n)
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
     &          + source(i,k,j)
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo

!       2nd: add in portion of vertical diffusion handled implicitly

        call ivdift (joff, js, je, istrt, iend, n, twodt)

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, js, je, istrt, iend, isst, itemp)
      call asbct (joff, js, je, istrt, iend, isss, isalt)
      call asbct (joff, js, je, istrt, iend, issdic, idic)
      call asbct (joff, js, je, istrt, iend, issdic13, idic13)
      call asbct (joff, js, je, istrt, iend, issalk, ialk)
      call asbct (joff, js, je, istrt, iend, isso2, io2)
      call asbct (joff, js, je, istrt, iend, isspo4, ipo4)
      call asbct (joff, js, je, istrt, iend, issphyt, iphyt)
      call asbct (joff, js, je, istrt, iend, isszoop, izoop)
      call asbct (joff, js, je, istrt, iend, issdetr, idetr)
      call asbct (joff, js, je, istrt, iend, isscocc, icocc)
      call asbct (joff, js, je, istrt, iend, isscaco3, icaco3)
      call asbct (joff, js, je, istrt, iend, issdop, idop)
      call asbct (joff, js, je, istrt, iend, issno3, ino3)
      call asbct (joff, js, je, istrt, iend, issdon, idon)
      call asbct (joff, js, je, istrt, iend, issdiaz, idiaz)
      call asbct (joff, js, je, istrt, iend, issdin15, idin15)
      call asbct (joff, js, je, istrt, iend, issdon15, idon15)
      call asbct (joff, js, je, istrt, iend, issphytn15, iphytn15)
      call asbct (joff, js, je, istrt, iend, isscoccn15, icoccn15)
      call asbct (joff, js, je, istrt, iend, isszoopn15, izoopn15)
      call asbct (joff, js, je, istrt, iend, issdetrn15, idetrn15)
      call asbct (joff, js, je, istrt, iend, issdiazn15, idiazn15)
      call asbct (joff, js, je, istrt, iend, issdfe, idfe)
      call asbct (joff, js, je, istrt, iend, issdetrfe, idetrfe)
      call asbct (joff, js, je, istrt, iend, issdoc13, idoc13)
      call asbct (joff, js, je, istrt, iend, issphytc13, iphytc13)
      call asbct (joff, js, je, istrt, iend, isscoccc13, icoccc13)
      call asbct (joff, js, je, istrt, iend, isscaco3c13, icaco3c13)
      call asbct (joff, js, je, istrt, iend, isszoopc13, izoopc13)
      call asbct (joff, js, je, istrt, iend, issdetrc13, idetrc13)
      call asbct (joff, js, je, istrt, iend, issdiazc13, idiazc13)

      return
      end

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, is, ie
      integer mask, m

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, dtdx, dtdy, dtdz
      real r2dt, cosdyt, fx, darea, boxar, rtwodt, sumdx, delx
      real sumdxr, dxdy, dxdydz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "ctavg.h"
      include "npzd.h"
      include "diag.h"
      include "diaga.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

      real temp1(imt,km), temp2(imt,km), temp3(imt,km)
      real twodt(km)

      include "isopyc.h"
      include "fdift.h"

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

      implicit none

      integer idiag, j, js, je, k, i, joff, iocv, jrow, is, ie

      real rdt, reltim, period

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "diaga.h"
      include "iounit.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "timeavgs.h"

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie

      real rts

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

      subroutine ivdift (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real rc2dt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"

      real twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif

      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)

!     compute residual implicit vertical mixing

      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif

      return
      end

      subroutine swflux0 (joff, js, je, is, ie, source)

      return
      end
!

