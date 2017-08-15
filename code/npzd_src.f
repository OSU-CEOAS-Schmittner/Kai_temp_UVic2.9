! source file: /data/home/kai/dev/UVic2.9/updates/npzd_src.F

      subroutine mobi_init

!     initialize MOBI parameters
      implicit none
      include "size.h"
      include "npzd.h"
      include "calendar.h"
      include "coord.h"
      include "stdunits.h"
      include "scalar.h"

      integer k, ioun, imobi
      real alpha, par
      real alpha_C, sumzpref

!     IO for fe_hydr
      include "param.h"
      include "pconst.h"
      integer ib(10), ic(10), iou, j, i
      character (120) :: fname, new_file_name
      logical exists
      real c1e3
      real tmpijkm(1:100,1:100,1:19)
! end IO for fe_hydr
      namelist /npzd/   alpha, kw, kc, abio_P, bbio, cbio, k1n, nup
     &,                 gamma1, gbio, epsbio, nuz, nud0
     &,                 wd0, par, dtnpzd, redctn, ki, redptn, capr
     &,                 dcaco3, nupt0, redotn, jdiar, kzoo, geZ
     &,                 zprefP, zprefZ, zprefDet, zprefDiaz, kfe, kfe_D
     &,                 kfemin, kfemax, kfemin_C, kfemax_C, kfeorg
     &,                 nupt0_D, sgbdfac, diazntp, nup_D, dfr, dfrt
     &,                 nudon0, nudop0, eps_assim, eps_excr, eps_nfix
     &,                 eps_wcdeni, eps_bdeni0, eps_recy, hdop
     &,                 mw, mwz, mw_c
     &,                 abio_C, k1n_C, nuct0,nuc, alpha_C
     &,                 zprefC, bapr
     &,                 wc0, kcal, dissk0, kc_c
!     set defaults for namelist npzd

      alpha = 0.16  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      kw = 0.04  ! Light attenuation due to water [1/m]
      kc = 0.047  ! Light atten. by phytoplankton [1/(m*mmol/m^3)]
      ki = 5.0  ! Light attenuation through sea ice & snow
      abio_P = 0.6  ! a; Maximum growth rate parameter [1/day]
      bbio = 1.066  ! b
      cbio = 1.0  ! c [1/deg_C]
      k1n = 0.7  ! Half saturation constant for N uptake [mmol/m^3]
      nup = 0.03  ! Specific mortality rate (Phytoplankton) [1/day]
      nup_D = 0.0001 !  Specific mortaility rate (Diazotrophs [1/day]
      nupt0 = 0.015  ! Fast-recycling  mortality rate (Phytoplankton) [1/day]
      nupt0_D = 0.001 ! Fast-recyling mortality rate (Diazotrophs) [1/day]
      gamma1 = 0.70  ! gama1; Assimilation efficiency (zpk)
      gbio = 0.38  ! Maximum grazing rate [1/day]
      epsbio = 1.6  ! Prey capture rate [(mmol/m^3)^(-2)day^(-1)]
      nuz = 0.06  ! Quadratic mortality (zpk) [(mmol/m^3)^(-2)day^(-1)]
      nud0 = 0.07  ! detritus remineralization rate [1/day]
      nudon0 = 2.33e-5 ! DON remineralization rate [1/day]
      nudop0 = 7.e-5  ! DOP remineralization rate [1/day]
      wd0 = 16.0  ! Sinking speed of detritus at surface [m/day]
      mwz = 100000. !Depth where sinking below remains constant (cm)
      mw = 0.02 !Sinking speed increase with depth (s-1)
      mw_c = 0.06 !Calcite sinking speed increase with depth (s-1)
      par = 0.43 ! fraction of photosythetically active radiation
      dtnpzd = dtts/4.  ! time step of biology [s]
      redctn = 7.1  ! C/N Redfield ratio (Schneider et al., 2003, GBC)
      redptn = 1./16.  ! P/N Redfield ratio
      capr = 0.022  ! carbonate to carbon production ratio
      dcaco3 = 650000.0  ! remineralisation depth of calcite [cm]
      redotn = 10.6  ! O2/N Redfield ratio (Anderson & Sarmiento, 1994, GBC)
      jdiar = 0.08  ! factor reducing the growth rate of diazotrophs
      kzoo = 0.15  ! half sat. constant for Z grazing [mmol N m^3]
      geZ = 0.6   ! Zooplankton growth efficiency
      sgbdfac = 1.0 ! sub-grid benthic denitrification rate factor
      diazntp = 28. !diazotroph N:P ratio
      dfr = 0.08   ! phyt mortality refractory/semi-labile DOM fraction
      dfrt = 0.01  ! phyt fast-recy refracotyr/semi-labile DOM fraction
      hdop = 0.4   ! DOP growth rate handicap
      alpha_C = 0.06  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      abio_C = 0.52  ! a; Maximum growth rate parameter coccolithophore [1/day]
      k1n_C = 0.4  ! Half sat const for N uptake (coccolithophores) [mmol/m^3]
      nuc = 0.03  ! Specific mortality rate (coccolithophores) [1/day]
      nuct0 = 0.015  ! Specific mortality rate (coccolithophores) [1/day]
      bapr = 0.05 ! detritus to carbonate ratio [mg POC/mg PIC]
      kc_c = 0.047  ! Light atten. by calcite [1/(m*mmol/m^3)]
      wc0 = 35. ! constant calcite sinking speed [m/day]
!      Gehlen et al 2007 use a surface sink speed of 50 m/day inc w/ depth
      kcal = 100.  !half sat for PIC dissolution from Planktom 10
      dissk0 = 0.013 !initial dissolution rate parameter [1/day]

      zprefP = 0.225 ! Zooplankton preference for P
      zprefC = 0.225 ! Zooplankton preference for C
      zprefZ = 0.225 ! Zooplankton preference for other Z
      zprefDet = 0.225 ! Zooplankton preference for Detritus
      zprefDiaz = 0.1 ! Zooplankton preference for diazotrophs
      eps_assim = 6.
      eps_excr = 4.
      eps_nfix = 1.
      eps_wcdeni = 25.
      eps_bdeni0 = 6.
      eps_recy = 1.
      kfemin = 0.04e-3    ! Minimum half saturation constant for Fe limitation [mmol Fe / m-3]
      kfemax = 0.2e-3    ! Maximum half saturation constant for Fe limitation [mmol Fe / m-3]
      pmax = 0.15    ! Phytoplankton biomass above which kfe increases [mmol N / m-3]
      kfe_D = 0.1e-3 ! Half saturation constant for Diaz Fe limitation [mmol Fe / m-3]
c     juan: Levin's value was 0.1e-3
      kfemin_C = 0.04e-3    ! Minimum half saturation constant for Fe limitation [mmol Fe / m-3]
      kfemax_C = 0.4e-3    ! Maximum half saturation constant for Fe limitation [mmol Fe / m-3]
      pmax_C = 0.15    ! Phytoplankton biomass above which kfe increases [mmol N / m-3]
      kfeleq = 10.**5.5   !  Fe-ligand stability constant [m^3/(mmol ligand)]
      lig = 1.0e-3      !  Ligand concentration  [mmol/m^3]
      thetamaxhi = 0.04    !  Maximum Chl:C ratio, abundant iron [gChl/(gC)]
      thetamaxlo = 0.01    !  Maximum Chl:C ratio, extreme iron limitation [gChl/(gC)]
      alphamax = 73.6e-6*86400  !   Maximum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
      alphamin = 18.4e-6*86400   !   Minimum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
      mc = 12.011           ! Molar mass of carbon [g/mol]
      fetopsed = 0.01  !  Fe:P for sedimentary iron source [molFe/molP]
      o2min = 5.      !  Minimum O2 concentration for aerobic respiration [mmolO_2/m^3]
      kfeorg = 0.45*(1./86400.)     !  Organic-matter dependent scavenging rate [(m^3/(gC s))^0.58]
      rfeton = 10.e-6*6.625        ! Uptake ratio of iron to nitrogen [mol Fe/mol N] = 10 micromol Fe /mol C
      kfecol = 0.005/86400.     ! Colloidal production and precipitation rate [s^-1]
!     read namelist
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, npzd, end=108)
108   continue
      write (stdout, npzd)
      call relunit (ioun)

!     convert units of NPZD parameters to MOM units
      redctn = redctn*1.e-3
      redotn = redotn*1.e-3
      redotp = redotn/redptn
      redctp = redctn/redptn
      redotc = redotn/redctn
      redntp = 1./redptn
      redntc = 1./redctn
      diazptn = 1./diazntp
      k1p_P   = k1n*redptn
      k1p_C  = k1n_C*redptn
      alpha_C = alpha_C/daylen
      abio_C = abio_C/daylen
      nuc = nuc/daylen
      nuct0 = nuct0/daylen
      wc0 = wc0*1.e2
      kc_c = kc_c*1.e-2
      dissk0 = dissk0/daylen
      kw = kw*1.e-2
      kc = kc*1.e-2
      ki = ki*1.e-2
      wd0 = wd0*1.e2
      alpha = alpha/daylen
      abio_P = abio_P/daylen
      nup = nup/daylen
      nup_D = nup_D/daylen
      nupt0 = nupt0/daylen
      nupt0_D = nupt0_D/daylen
      gbio = gbio/daylen
      epsbio = epsbio/daylen
      nuz = nuz/daylen
      nud0 = nud0/daylen
      nudop0 = nudop0/daylen
      nudon0 = nudon0/daylen
      tap = 2.*alpha*par
      tap_C = 2.*alpha_C*par !Total Attenuation by Phyt?

!     calculate sinking speed of detritus divided by grid width
      do k=1,km
!     linear increase wd0-200m with depth
         if (zt(k) .lt. mwz) then
            wd(k) = (wd0+mw*zt(k))/daylen/dzt(k) ! [s-1]
            wc(k) = (wc0+mw_c*zt(k))/daylen/dzt(k)    ! [s-1]
         else
            wd(k) = (wd0+mw*mwz)/daylen/dzt(k) ! [s-1]
            wc(k) = (wc0+mw_c*mwz)/daylen/dzt(k)    ! [s-1]
         endif
         rkwz(k) = 1./(kw*dzt(k))
      enddo
      ztt(1)=0.0
      do k=1,km
         ztt(k+1)=(-1)*zw(k)
      enddo

!     check grazing preferences = 1 for N case
      sumzpref = zprefP + zprefDet + zprefZ + zprefDiaz + zprefC
      IF (sumzpref.ne.1.) THEN
         print*,"mobi_init: adjust zprefs from"
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefDiaz, zprefC
         zprefP = zprefP/sumzpref
         zprefC = zprefC/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         zprefDiaz = zprefDiaz/sumzpref
         print*,"to"
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefDiaz, zprefC
      END IF
! read the hydrothermal Fe input from O_fe_hydr.nc
      fe_hydr(:,:,:) = 0.
      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = 19

      fname = new_file_name ("O_fe_hydr.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         c1e3 = 1000
         call openfile (trim(fname), iou)
         call getvara ('O_fe_hydr', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1e3, c0)
         fe_hydr(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)

            do k=1,19
               do j=1,jmt
                  fe_hydr(1,j,k) = fe_hydr(imtm1,j,k)
                  fe_hydr(imt,j,k) = fe_hydr(2,j,k)
               enddo
               do i=1,imt
                  fe_hydr(i,1,k) = fe_hydr(i,2,k)
                  fe_hydr(i,jmt,k) = fe_hydr(2,j,k)
               enddo
            enddo
      else
         print*,"Warning => Cannot find", trim(fname)
      endif

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = 1./dzt(1)
      do k=2,km
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
      imobi=1
      call setimobi('imobipo4',imobipo4,imobi)
      call setimobi('imobiphyt',imobiphyt,imobi)
      call setimobi('imobizoop',imobizoop,imobi)
      call setimobi('imobidetr',imobidetr,imobi)
      call setimobi('imobidic',imobidic,imobi)
      call setimobi('imobidic13',imobidic13,imobi)
      call setimobi('imobiphytc13',imobiphytc13,imobi)
      call setimobi('imobizoopc13',imobizoopc13,imobi)
      call setimobi('imobidetrc13',imobidetrc13,imobi)
      call setimobi('imobidoc13',imobidoc13,imobi)
      call setimobi('imobidiazc13',imobidiazc13,imobi)
      call setimobi('imobicoccc13',imobicoccc13,imobi)
      call setimobi('imobicaco3c13',imobicaco3c13,imobi)
      call setimobi('imobidop',imobidop,imobi)
      call setimobi('imobino3',imobino3,imobi)
      call setimobi('imobidon',imobidon,imobi)
      call setimobi('imobidiaz',imobidiaz,imobi)
      call setimobi('imobidin15',imobidin15,imobi)
      call setimobi('imobidon15',imobidon15,imobi)
      call setimobi('imobiphytn15',imobiphytn15,imobi)
      call setimobi('imobizoopn15',imobizoopn15,imobi)
      call setimobi('imobidetrn15',imobidetrn15,imobi)
      call setimobi('imobidiazn15',imobidiazn15,imobi)
      call setimobi('imobicoccn15',imobicoccn15,imobi)
      call setimobi('imobicocc',imobicocc,imobi)
      call setimobi('imobicaco3',imobicaco3,imobi)
      call setimobi('imobidfe',imobidfe,imobi)
      call setimobi('imobidetrfe',imobidetrfe,imobi)
      if (imobi-1 .eq. ntnpzd) then
         print*,'total number of MOBI tracers =', imobi-1
      else
         print*,'error in mobi_init:'
         print*,'imobi-1=',imobi-1,' not equal to ntnpzd', ntnpzd
         stop
      endif
      return
      end
!     END mobi_init

      subroutine setimobi (name, index, inc)
      implicit none
      character(*) :: name
      integer index, inc
      index = inc
      inc = index + 1
      print*,name,' = ',index
      return
      end

      subroutine mobi_driver(
     &                 kmx, twodt, rctheta, dayfrac, swr, tnpzd
     &,                t_in, po4_in
     &,                o2_in
     &,                s_in, dic_in, alk_in, co2_in
     &,                dic13_in
     &,                sgb_in
     &,                no3_in
     &,                din15_in
     &,                src)

!     main driver for Model of Ocean Biogeochemistry and Isotopes (MOBI)
!     input: kmx, twodt, rctheta, dayfrac, swr, tnpzd, t_in,
!            po4_in, (o2_in, s_in, dic_in, alk_in, co2_in, dic13_in, no3_in,
!            sgb_in, din15_in)
!     output: src, (sedcorgflx, sedcalflx)
!             if 1 is switch on additional output (rnpp etc.)
!                is in npzd.h
!     copied from tracer.F (Andreas Schmittner, Oct 2, 2013)
!     modified to a vertical column model
      implicit none
      integer k,kmx
      include "size.h"
      include "coord.h"
      include "grdvar.h"
      include "mw.h"
      include "npzd.h"

      real snpzd(ntnpzd), tnpzd(km,ntnpzd), src(km,nsrc)
      real dic_npzd_sms(km), twodt, rctheta, gl, impo, expo
      real npp, remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real graz_Det, graz_Z, phin, dz, prca, dprca, nud, bct, bctz
      real t_in(km), po4_in(km), sedrr, nudop, nudon, npp_dop
      real impofe, expofe, remife, nfix(km)
      real fesed, bfe, sgb, o2_in(km)
      real fo2, so2, o2flag
      real calpro
      real biococc,jmax_C,u_C,avej_C,gd_C,npp_C, morp_C,morpt_C
      real nuct, g_C, graz_C, ing_C, felimit_C, gl_C, npp_C_dop
      real dissl, caco3in, sil_in, atmpres1, pHlo1, pHhi1, pH1, p_in
      real co2star1, dco2star1, pCO21, dpco21, CO31, omegaca, omegaar
      real c_in, wwc,expocaco3, impocaco3, dissk1, calatt
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
      real rc13impo, rc13expo, ac13_DIC_aq
      real ac13_aq_POC, ac13b, prca13, rdic13, rtdic13
      real pH, co2star, dco2star, pCO2, dpco2, CO3
      real Omega_c, Omega_a, pCO2_old, pCO2_new
      real atmpres, depth
      real rcaco3c13impo, rcaco3c13expo, rtcaco3c13
      real sgb_in(km)
      real no3_in(km), lno3, sg_bdeni
      real morp_D, lntp, npp_D_dop
      real npp_D, graz_D, morpt_D, no3flag, wcdeni, bdeni(km)
      real din15_in(km), din15flag, eps_bdeni
      real rn15impo, rn15expo, uno3, rno3, bwcdeni, bbdeni

      expo = 0.0
      impo = 0.0
      phin = 0.0                ! integrated phytoplankton
      prca = 0.0                ! integrated production of calcite

      sedrr = 0.0
      rsedrr = 0.0
      rprca = 0.0
      rn15impo = 0.0
      rn15expo = 0.0
      rc13impo = 0.0
      rc13expo = 0.0
      prca13 = 0.0
      rcaco3c13impo = 0.0
      rcaco3c13expo = 0.0

      expofe = 0.0
      impofe = 0.0
      calpro = 0.0
      caco3in = 0.0
      calatt = 0.0
      dissl = 0.0
      impocaco3 = 0.0
      expocaco3 = 0.0
      dissk1 = 0.0

!1111 start of main k-loop
      do k=1,kmx
         rn15impo = rn15expo
         atmpres = 1.0          !atm
         depth = zt(k)/100.     !m

         call co2calc_SWS (t_in(k), s_in(k), dic_in(k), alk_in(k)
     &,                    co2_in, atmpres, depth, pH, co2star, dco2star
     &,                    pCO2, dpco2
     &,                    CO3, Omega_c, Omega_a)
!        c13 biological fractionation
         ac13_DIC_aq = -1.0512994e-4*t_in(k)+1.011765
!        Popp et al. (1989) Am. J. Sci.
         ac13_aq_POC = -0.017*log10(min(max(co2star*1000.,2.),74.))
     &                 +1.0034
         ac13b = ac13_aq_POC/ac13_DIC_aq
         rc13impo = rc13expo
         rcaco3c13impo = rcaco3c13expo
c         dissk0 =(1-(CO3-CO3/Omega_c))/(kcal+(CO3-CO3/Omega_C))
c         dissk1 = min(1.,max(0.,dissk0)) ! Andreas
         dissk1 = dissk0*(1-Omega_c) ! Andreas 11 d^-1 from Gehlen 2006
         swr = swr*exp(-kc*phin
     &             - kc_c*caco3in
     &              )

         phin = max(tnpzd(k,imobiphyt), trcmin)*dzt(k)
     &        + max(tnpzd(k,imobidiaz), trcmin)*dzt(k)
     &        + max(tnpzd(k,imobicocc), trcmin)*dzt(k)
         caco3in = caco3in + tnpzd(k,imobicaco3)*dzt(k)
         impocaco3 = expocaco3*dztr(k)
         gl = swr*exp(ztt(k)*rctheta)
         impo = expo*dztr(k)
         impofe = expofe*dztr(k)

         bct = bbio**(cbio*t_in(k))
         bctz = (0.5*(tanh(o2_in(k) - 8.)+1))
     &        *bbio**(cbio*min(t_in(k),20.))
!        decrease remineralisation rate in oxygen minimum zone
         nud = nud0*(0.65+0.35*tanh(o2_in(k) - 3.0))
         nudon = nudon0
         nudop = nudop0
!-----------------------------------------------------------------------
!        call the npzd ecosystem dynamics model
!-----------------------------------------------------------------------

         call mobi_src (tnpzd(k,:), gl, bct, impo, dzt(k)
     &,                 dayfrac, wd(k), rkwz(k), nud
     &,                 impocaco3,wc(k),dissk1, expocaco3
     &,                 nudop, nudon, nfix(k)
     &,                 snpzd, expo
     &,                 calpro
     &,                 npp, morpt, remi, excr, graz, morp, morz
     &,                 graz_Det, graz_Z
     &,                 npp_D, npp_D_dop, graz_D, morp_D, morpt_D
     &,                 npp_dop
     &,                 npp_C_dop
     &,                 npp_C, morpt_C, graz_C
     &,                 morp_C, dissl, calatt
     &,                 bctz
     &,                 rn15impo, rn15expo
     &,                 rc13impo, rc13expo, ac13b
     &,                 rcaco3c13impo, rcaco3c13expo
     &,                    expofe, impofe, remife
c     juan: use o2_in(k) instead of t(i,k,j,io2,taum1)
     &,                    o2_in(k)
     &                      )
! These are source/sink terms
         snpzd = snpzd*rdtts
         expofe = expofe*rnbio
         expocaco3 = expocaco3*rnbio
         expo = expo*rnbio
         rn15expo = rn15expo*rnbio
         rc13expo = rc13expo*rnbio
         rcaco3c13expo = rcaco3c13expo*rnbio
         rcalpro(k) = calpro*rnbio
         rcalatt(k) = calatt*rnbio
         rdissl(k) = dissl*rdtts
         rexpo(k) = expo
         rgraz(k) = graz*rnbio
         rgraz_Det(k) = graz_Det*rnbio
         rgraz_Z(k) = graz_Z*rnbio
         rmorp(k) = morp*rnbio
         rmorz(k) = morz*rnbio
         rnpp(k) = npp*rnbio
         rmorpt(k) = morpt*rnbio
         rremi(k) = remi*rnbio
         rexcr(k) = excr*rnbio
         rnpp_dop(k) = npp_dop*rnbio
         rnpp_D(k) = npp_D*rnbio
         rnpp_D_dop(k) = npp_D_dop*rnbio
         rgraz_D(k) = graz_D*rnbio
         rmorp_D(k) = morp_D*rnbio
         rmorpt_D(k) = morpt_D*rnbio
         rnfix(k) = nfix(k)*rnbio
         rnpp_C_dop(k) = npp_C_dop*rnbio
         rexpocaco3(k) = expocaco3
         rnpp_C(k) = npp_C*rnbio
         rgraz_C(k) = graz_C*rnbio
         rmorp_C(k) = morp_C*rnbio
         rmorpt_C(k) = morpt_C*rnbio
         rexpofe(k) = expofe
         rremife(k) = remife*rnbio
!----------------------------------------------------------------------
!             calculate calcite at the bottom and dissolve
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

         sedrr = sgb_in(k)*expo*dzt(k)
!-----------------------------------------------------------------------
!        benthic denitrification model of Bohlen et al., 2012, GBC (eq. 10)
!        NO3 is removed out of bottom water nitrate.
!        See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!        limit denitrification as nitrate approaches 0 uM
         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
         lno3 = 0.5*tanh(no3_in(k)*10 - 5.0)

         sg_bdeni = (0.06 + 0.19*0.99**(max(o2_in(k),trcmin)
     &                 - max(no3_in(k),trcmin)))
     &              *max(expo*sgb_in(k),trcmin)*redctn*1.e3
         sg_bdeni = min(sg_bdeni, sgb_in(k)*expo)
         sg_bdeni = max(sg_bdeni, 0.)

c Andreas is this line correct? Why 2000 m?
c         if (zt(k) .le. 200000.) sg_bdeni = sgbdfac*sg_bdeni

         sg_bdeni = sg_bdeni*(0.5 + lno3)*no3flag*din15flag
         bdeni(k) = sg_bdeni
         snpzd(imobino3) = snpzd(imobino3) + sgb_in(k)*expo
     &                                           - sg_bdeni
         rbdeni(k) = bdeni(k)
         rno3 = max(din15_in(k), trcmin*rn15std/(1+rn15std))
     &            /max(no3_in(k)-din15_in(k),trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         eps_bdeni = eps_bdeni0*exp(-2.5e-6*(zt(k)))
         bbdeni = rno3 - eps_bdeni*rno3/1000.

         snpzd(imobidin15) = snpzd(imobidin15)
     &                     + rn15expo*sgb_in(k)*expo
     &                     - bbdeni/(1+bbdeni)*sg_bdeni
         fesed=fetopsed*bct*redptn*expo*sgb_in(k)
         bfe=fesed
         snpzd(imobidfe) = snpzd(imobidfe) + fesed
         o2flag = 0.5 - sign(0.5, o2_in(k)-o2min)
         snpzd(imobidfe) = snpzd(imobidfe) + expofe*sgb_in(k)*o2flag
         bfe = bfe + expofe*sgb_in(k)*o2flag
         expofe = expofe - sgb_in(k)*expofe*o2flag
         rremife(k) = rremife(k) + bfe
         snpzd(imobipo4) = snpzd(imobipo4) + sgb_in(k)*expo*redptn
         snpzd(imobidic) = snpzd(imobidic) + sgb_in(k)*expo*redctn
         snpzd(imobidic13) = snpzd(imobidic13)
     &                     + rc13expo*sgb_in(k)*expo*redctn
         expo = expo - sgb_in(k)*expo
         rremi(k) = rremi(k) + sgb_in(k)*expo
         rsedrr = rsedrr + sedrr

!-----------------------------------------------------------------------
!             set source/sink terms
!-----------------------------------------------------------------------

         src(k,ispo4) = snpzd(imobipo4)
         src(k,isphyt) = snpzd(imobiphyt)
         src(k,iszoop) = snpzd(imobizoop)
         src(k,isdetr) = snpzd(imobidetr)
         src(k,isdic) = snpzd(imobidic)
         src(k,isdop) = snpzd(imobidop)
         src(k,isno3) = snpzd(imobino3)
         src(k,isdon) = snpzd(imobidon)
         src(k,isdiaz) = snpzd(imobidiaz)
         src(k,isdin15) = snpzd(imobidin15)
         src(k,isdon15) = snpzd(imobidon15)
         src(k,isphytn15) = snpzd(imobiphytn15)
         src(k,iszoopn15) = snpzd(imobizoopn15)
         src(k,isdetrn15) = snpzd(imobidetrn15)
         src(k,isdiazn15) = snpzd(imobidiazn15)
         src(k,iscoccn15) = snpzd(imobicoccn15)
         src(k,isdic13) = snpzd(imobidic13)
         src(k,isphytc13) = snpzd(imobiphytc13)
         src(k,iszoopc13) = snpzd(imobizoopc13)
         src(k,isdetrc13) = snpzd(imobidetrc13)
         src(k,isdoc13) = snpzd(imobidoc13)
         src(k,isdiazc13) = snpzd(imobidiazc13)
         src(k,iscoccc13) = snpzd(imobicoccc13)
         src(k,iscaco3c13) = snpzd(imobicaco3c13)
         src(k,iscocc) = snpzd(imobicocc)
         src(k,iscaco3) = snpzd(imobicaco3)
         src(k,isdfe) = snpzd(imobidfe)
         src(k,isdetrfe) = snpzd(imobidetrfe)
!        These are organic sources and sinks of DIC (i.e. remin - pp)
!        all are based on po4 uptake and remineralization
         dic_npzd_sms(k) = snpzd(imobidic)
!        production of calcite
         dprca = rcalpro(k)*1e-3
         prca = prca + dprca*dzt(k)

         rtdic13 = max(dic13_in(k),trcmin*rc13std/(1+rc13std))
     &             / max(dic_in(k),trcmin)
         rtdic13 = min(rtdic13, 2.*rc13std/(1+rc13std))
         rtdic13 = max(rtdic13, 0.5*rc13std/(1+rc13std))

         prca13 = prca13 + dprca*dzt(k)*rtdic13
         rtcaco3c13 = max(tnpzd(k,imobicaco3c13)
     &                     ,trcmin*rc13std/(1+rc13std))
     &             / max(tnpzd(k,imobicaco3c13),trcmin)
         rtcaco3c13 = min(rtcaco3c13, 2.*rc13std/(1+rc13std))
         rtcaco3c13 = max(rtcaco3c13, 0.5*rc13std/(1+rc13std))
         src(k,isalk) = -snpzd(imobidic)*redntc*1.e-3
!        calculate total export to get total import for next layer
         expo = expo*dzt(k)
         expofe = expofe*dzt(k)
         expocaco3 = expocaco3*dzt(k)
      enddo
!1111 end of main k-loop
      rprca = prca

!2222 start of second k-loop
      do k=1,kmx
!        limit oxygen consumption below concentrations of
!        5umol/kg as recommended in OCMIP
         fo2 = 0.5*tanh(o2_in(k) - 2.5)
!        sink of oxygen
         so2 = dic_npzd_sms(k)*redotc
! O2 is needed to generate the equivalent of NO3 from N2 during N2 fixation
! 0.5 H2O + 0.5 N2+1.25O2 -> HNO3
! note that so2 is -dO2/dt
     &          + nfix(k)*rnbio*1.25e-3
!        add denitrification as source term for NO3
         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
         lno3 = 0.5*tanh(no3_in(k) - 2.5)
         lntp = 0.5*tanh(no3_in(k)/(redntp*po4_in(k))*100 - 60.)
!        800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
         wcdeni = 800.*no3flag*so2*(0.5 - fo2)*(0.5 + lno3)
     &           *din15flag
         wcdeni = max(wcdeni, 0.)
         src(k,isno3) = src(k,isno3) - wcdeni
!        calculate isotope effect of water column denitrification
         uno3 = wcdeni*twodt/no3_in(k)
         uno3 = min(uno3, 0.999)
         uno3 = max(uno3, trcmin)
         rno3 = max(din15_in(k),trcmin*rn15std/(1+rn15std))
     &             / max(no3_in(k)-din15_in(k),
     &                   trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         bwcdeni = rno3 + eps_wcdeni*(1-uno3)/uno3*log(1-uno3)*rno3/1000.

         src(k,isdin15) = src(k,isdin15) - (bwcdeni/(1+bwcdeni))*wcdeni
! Correct the ALK stoichiometry to account for N2 fixation
         src(k,isalk) = src(k,isalk) + wcdeni*1.e-3
!bdeni
         src(k,isalk) = src(k,isalk) + bdeni(k)*1.e-3
! Now add account for N2 fixation (ALK production is tied to PO4 change and
! thus in the case of N2 fixation is not correct).
         src(k,isalk) = src(k,isalk) - nfix(k)*rnbio*1.e-3
         rwcdeni(k) = wcdeni
         src(k,iso2) = -so2*(0.5 + fo2)
      enddo
!2222 end of second k-loop

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
!3333 start of third k-loop
      do k=1,kmx-1
         src(k,isdic) = src(k,isdic)
     &        + rdissl(k)*1.e-3
     &        - rcalpro(k)*1.e-3
     &        - rcalatt(k)*1.e-3
         src(k,isdic13) = src(k,isdic13)
     &        + rdissl(k)*1.e-3*rtcaco3c13 ! AS
     &        - rcalpro(k)*1.e-3*rtdic13
     &        - rcalatt(k)*1.e-3*rtdic13
         src(k,isalk) = src(k,isalk)
     &        + 2.*rdissl(k)*1.e-3
     &        - 2.*rcalpro(k)*1.e-3
     &        - 2.*rcalatt(k)*1.e-3
      enddo
!3333 end of third k-loop
      src(kmx,isdic) = src(kmx,isdic)
     &     + rdissl(kmx)*1.e-3
     &     - rcalpro(kmx)*1.e-3
     &     - rcalatt(kmx)*1.e-3
     &     + rexpocaco3(kmx)*1.e-3
      src(kmx,isdic13) = src(kmx,isdic13)
     &                 + rdissl(kmx)*1.e-3*rtcaco3c13
     &                 - rcalpro(kmx)*1.e-3*rtdic13
     &                 - rcalatt(kmx)*1.e-3*rtdic13
     &                 + rexpocaco3(kmx)*1.e-3*rtdic13
      src(kmx,isalk) = src(kmx,isalk)
     &               + 2.*rdissl(kmx)*1.e-3
     &               - 2.*rcalpro(kmx)*1.e-3
     &               - 2.*rcalatt(kmx)*1.e-3
     &               + 2.*rexpocaco3(kmx)*1.e-3
      return
      end
!     END mobi_driver

      subroutine mobi_src (bioin, gl, bct, impo, dzt
     &,                    dayfrac, wwd, rkw, nud
     &,                    impocaco3, wwc, dissk1, expocaco3out
     &,                    nudop, nudon, nfixout
     &,                    bioout, expoout
     &,                    calproout
     &,                    nppout, morptout, remiout, excrout, grazout
     &,                    morpout, morzout, graz_Det_out, graz_Zout
     &,                    npp_Dout, npp_D_dopout, graz_Dout, morp_Dout
     &,                    morpt_Dout, npp_dopout
     &,                    npp_C_dopout
     &,                    npp_Cout, morpt_Cout,graz_Cout
     &,                    morp_Cout, disslout, calattout
     &,                    bctz
     &,                    rn15impo, rn15expoout
     &,                    rc13impo, rc13expoout, ac13b
     &,                    rcaco3c13impo, rcaco3c13expoout
     &,                    expofeout, impofe, remifeout
     &,                    o2
     &                     )

!=======================================================================
!     computes source terms of the NPZD model
!     initial version of code adapted from Xavier Giraud:
!     Giraud et al. 2000, J Mar Res, 58, 609-630
!     original model reference:
!     Oeschlies and Garcon 1999, Global Biogeochem. Cycles 13, 135-160
!     Schmittner et al. 2005,  Global Biogeochem. Cycles 19, GB3004,
!     doi:10.1029/2004GB002283.
!     Schmittner et al. 2008, Global Biogeochem. Cycles 22, GB1013
!
!     This version was modified by David Keller and corrects the zooplankton
!     grazing formulation.  Note that zooplankton are now allowed to graze
!     on themselves and detritus, in addition to phyt. and diazotrophs.
!     The calculation of light has also been corrected.
!     Keller et al. 2012, Geosci. Model Dev., 5, 1195-1220
!
!     The nitrogen isotope model (nitrogen_15) has been developed by
!     Chris Somes and is described in
!     Somes et al. 2010, Global Biogeochem. Cycles 24, GB4019
!     Somes et al. 2010, Geophys. Res. Lett. 37, L23605 and
!     Somes et al. 2013, Biogeosc. 10, 5889-5910
!
!     The carbon isotope model (carbon_13) has been written by Andreas
!     Schmittner and is documented in
!     Schmittner et al. 2013 Biogeosc. 10, 5793-5816
!     Chris Somes modified the code converting from alpha to beta formulation
!
!     Karen Kvale has included prognostic equations for CaCO3, coccolithphores, and
!     ballast described in
!     Kvale et al. 2015 Atmos.-Ocean 53, doi:10.1080/07055900.2015.1049112
!
!     Levin Nickelsen has written the iron model:
!     Nickelsen et al. 2015, Geosc. Model Dev. 8, 1357-1381, doi:10.5194/gmd-8-1357-2015
!
!     Juan Muglia has modified the subgrid-scale scheme used for benthic denitrification
!     and iron release from sediments. He has also included geothermal sources of iron.
!
!     Note that nutrient (N) represents phosphate
!
!     input variables:
!
!       bioin(1:ntnpzd) [mmol m-3]
!
!       gl         = 2.*light at top of grid box
!       nbio       = number of time steps
!       dtbio        = time step [s]
!       bct        = bbio**(cbio*temperature)
!       impo       = import of detritus from above [mmol m-3]
!       dzt        = depth of grid box [cm]
!       dayfrac    = day length (fraction: 0 < dayfrac < 1)
!       wwd        = sinking speed of detritus/dzt
!       rkw        = reciprical of kw*dzt(k)
!       nud        = remineralisation rate of detritus [s-1]
!
!     output variables:
!
!       bioout     = change from bioin [mmol m-3]
!       nppout     = net primary production [mmol m-3]
!       grazout    = grazing [mmol m-3]
!       morpout    = quadratic mortality of phytoplankton [mmol m-3]
!       morptout   = specific mortality of phytoplankton [mmol m-3]
!       morzout    = mortality of zooplankton [mmol m-3]
!       remiout    = remineralisation [mmol m-3]
!       excrout    = excretion [mmol m-3]
!       expoout    = detrital export [mmol m-3]
!       npp_Dout   = NPP of diazotrophs
!       graz_Dout  = grazing of diazotrophs
!       morp_Dout  = mortality of diazotrophs
!       nfixout    = rate of N2 fixation
!       graz_Det_out = grazing of detritus
!       graz_Zout   = grazing on othe zooplankton
!       avej_out    = light-depend phyt. growth rate
!       avej_D_out  = light-depend Diaz growth rate
!       gmax_out    = temp-depend. zoo growth rate
!       no3P_out    = no3 depend. phyt growth rate
!       po4P_out    = po4 depend. phyt growth rate
!       po4_D_out   = po4 depend. Diaz growth rate
!
!       New grazing formulation variables and parameters
!
!       The following terms determine ingestion according to a
!       a Holling II curve (i.e. Michaelis Menten):
!
!       Ingestion = max_graz_rate * (Ft/(Ft + kzoo))
!
!       where Ft is the weighted measure of the total food available
!       and equals the sum of the different prey types times the
!       preference of Z for that type of prey
!
!       zprefP   = Z preference for P
!       zprefDiaz   = Z preference for Diaz
!       zprefDet = Z preference for detritus
!       zprefZ   = Z preference for other Z
!       kzoo = half saturation coefficienct for Z ingestion mmol N m-3
!       ing_P    = zooplankton ingestion of phytoplankon
!       ing_D    = zooplankton ingestion of diazotrophs
!       ing_Det  = zooplankton ingestion of detritus
!       ing_Z    = zooplankton ingestion of other zooplankton
!       thetaZ   = Michaelis-Menten denominator
!
!=======================================================================

      implicit none

      integer n

      real gl, f1, biopo4, biophyt, biozoop, biodetr, jmax, u_P, g_P
      real morp, morpt, morz, remi, excr, expo, impo, nppout, grazout
      real morpout, morptout, morzout, remiout, excrout, expoout
      real avej_out, avej_D_out, gmax_out, no3P_out, po4P_out, po4_D_out
      real dzt, po4flag, phytflag, zoopflag, detrflag, wwd, rkw, gd
      real nupt, nud, biono3, u_D,npp_D, npp_Dout, no3flag, biodiaz, bct
      real diazflag, g_D,graz_D, morpt_D, jmax_D, gd_D, avej_D, no3upt_D
      real morpt_Dout, graz_Dout, nfixout, biop2, u1, u2, phi1, phi2
      real avej, graz_Det_out, graz_Zout, thetaZ, ing_P, ing_D, dayfrac
      real ing_Det, ing_Z, g_Z, g_Det, graz_Z, graz_Det, gmax
      real no3P, po4P, po4_D, bctz
      real excrdiaz, biodop, biodon, dopflag, donflag, recy_don
      real recy_dop, npp, graz, biodin15, biodon15, biophytn15
      real biodetrn15, biodiazn15, bassim, fcassim, bexcr, fcexcr
      real bnfix, fcnfix, rn15impo, rn15expoout, dig, dig_P, dig_Z
      real dig_Det, dig_D, excr_P, excr_Z, excr_Det, excr_D, sf, sf_P
      real sf_Z, sf_Det, sf_D, nr_excr_D, uno3, rno3, rzoop
      real rtdin15, rtphytn15, rtzoopn15, rtdetrn15, rtdiazn15, rtdon15
      real din15flag, biozoopn15, morp_Dout
      real don15flag, phytn15flag, zoopn15flag, detrn15flag, diazn15flag
      real limP, limP_D, dopupt_D_flag, dopupt_D, morp_D, nupt_D
      real limP_dop, limP_po4
      real udon, rdon, brecy, fcrecy, biodic, dopupt_flag, dopupt
      real biodic13, biophytc13, biozoopc13, biodetrc13, biodiazc13
      real dic13flag, phytc13flag, zoopc13flag, detrc13flag, diazc13flag
      real biodoc13, doc13flag, biodoc, rcaco3c13impo, rcaco3c13expoout
      real rdic13, rtphytc13, rtzoopc13, rtdetrc13, rtdiazc13, ac13b
      real rc13impo, rc13expoout, dicflag, bc13npp, rtdic13, fcnpp
      real rtdoc13, nudop, nudon, npp_dopout, npp_D_dopout

      real expo_B, impo_B, remi_B, remiout_B, expoout_B, dflag_B
      real graz_Det_out_B, graz_Det_B, g_Det_B, biod_B, dig_Det_B
      real graz_Z_B, graz_C_B, morp_C_B, morz_B, excr_Det_B, sf_Det_B

      real biococc,jmax_C,u_C,avej_C,gd_C,npp_C,npp_Cout,morp_C,morpt_C
      real nuct,g_C,coccflag,graz_C,graz_Cout,morp_Cout,coccn15flag
      real morpt_Cout, ing_C, gl_C, limP_C, dig_C, excr_C, sf_C
      real dopupt_C_flag, dopupt_C, npp_C_dopout, biococcn15,rtcoccn15
      real biococcc13, coccc13flag, rtcoccc13, biocaco3c13, caco3c13flag
      real rtcaco3c13

      real caco3flag, wwc,expocaco3
      real expocaco3out,calpro,calproout
      real biocaco3, impocaco3, dissk1, dissl, disslout
      real calatt,calattout

      real biodfe, biodetrfe, dfeflag, detrfeflag, o2flag
      real expofe, impofe, feorgads, remife, thetamax, deffe, fepa
      real thetamax_D, deffe_D, thetachl, thetachl_D, chl, feprime
      real fecol, irrtop, kirr, aveirr, fesed, gl_O, gl_D, alpha_O
      real alpha_D, par, deffe_C, thetamax_C, thetamaxout_C, deffeout_C

      real expofeout, impofeout, feorgadsout, fecolout
      real remifeout, thetamaxout, deffeout, thetachlout
      real chlout, feprimeout, o2, chl_D_out

      real p1,p2,kfevar, p1_C, p2_C, kfevar_C

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "npzd.h"

      real bioin(ntnpzd), bioout(ntnpzd)

      p1=min(bioin(imobiphyt),pmax) !phytoplankton
      p2=max(0.0,bioin(imobiphyt)-pmax)
      kfevar=(kfemin*p1+kfemax*p2)/(p1+p2)
      deffe=bioin(imobidfe)/(kfevar+bioin(imobidfe))
      thetamax=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe
      alpha_O=alphamin+(alphamax-alphamin)*deffe
      gl_O=gl*thetamax*alpha_O
      p1_C=min(bioin(imobicocc),pmax_C) !calcifiers
      p2_C=max(0.0,bioin(imobicocc)-pmax_C)
      kfevar_C=(kfemin_C*p1_C+kfemax_C*p2_C)/(p1_C+p2_C)
      deffe_C=bioin(imobidfe)/(kfevar_C+bioin(imobidfe))
      thetamax_C=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe_C
      alpha_O=alphamin+(alphamax-alphamin)*deffe_C
      gl_C=gl*thetamax_C*alpha_O
      deffe_D=bioin(imobidfe)/(kfe_D+bioin(imobidfe))
      thetamax_D=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe_D
      alpha_D=alphamin+(alphamax-alphamin)*deffe_D
      gl_D=gl*thetamax_D*alpha_D

      bioout(:) = 0.0
      biopo4 = bioin(imobipo4)
      biophyt = bioin(imobiphyt)
      biozoop = bioin(imobizoop)
      biodetr = bioin(imobidetr)
      biodic = bioin(imobidic)
      biodop = bioin(imobidop)
      biono3 = bioin(imobino3)
      biodon = bioin(imobidon)
      biodiaz = bioin(imobidiaz)
      biodin15 = bioin(imobidin15)
      biodon15 = bioin(imobidon15)
      biophytn15 = bioin(imobiphytn15)
      biozoopn15 = bioin(imobizoopn15)
      biodetrn15 = bioin(imobidetrn15)
      biodiazn15 = bioin(imobidiazn15)
      biococcn15 = bioin(imobicoccn15)
      biodic13 = bioin(imobidic13)
      biophytc13 = bioin(imobiphytc13)
      biozoopc13 = bioin(imobizoopc13)
      biodetrc13 = bioin(imobidetrc13)
      biodoc13 = bioin(imobidoc13)
      biodiazc13 = bioin(imobidiazc13)
      biococcc13 = bioin(imobicoccc13)
      biocaco3c13 = bioin(imobicaco3c13)
      biococc = bioin(imobicocc)
      biocaco3 = bioin(imobicaco3)
      biodfe = bioin(imobidfe)
      biodetrfe = bioin(imobidetrfe)

!       flags prevent negative values by setting outgoing fluxes to
!     zero if tracers are lower than trcmin
      po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
      phytflag = 0.5 + sign(0.5,biophyt - trcmin)
      zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
      detrflag = 0.5 + sign(0.5,biodetr - trcmin)
!     set defaults
      no3flag = 1.
      dopflag = 1.
      donflag = 1.
      diazflag = 1.
      din15flag = 1.
      don15flag = 1.
      phytn15flag = 1.
      zoopn15flag = 1.
      detrn15flag = 1.
      diazn15flag = 1.
      dic13flag = 1.
      doc13flag = 1.
      phytc13flag = 1.
      zoopc13flag = 1.
      detrc13flag = 1.
      diazc13flag = 1.
      dfeflag = 1.
      coccflag = 1.
      detrfeflag = 1.
      coccn15flag = 1.
      coccc13flag = 1.
      caco3c13flag = 1.
      dopflag = 0.5 + sign(0.5,biodop - trcmin)
      no3flag = 0.5 + sign(0.5,biono3 - trcmin)
      donflag = 0.5 + sign(0.5,biodon - trcmin)
      diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
      din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
      don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
      phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
      coccn15flag = 0.5 + sign(0.5, biococcn15 - trcmin)
      zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
      detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
      diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
      dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
      phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
      coccc13flag = 0.5 + sign(0.5,biococcc13 - trcmin)
      caco3c13flag = 0.5 + sign(0.5,biocaco3c13 - trcmin)
      zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
      detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
      doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
      diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
      dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
      detrfeflag = 0.5 + sign(0.5,biodetrfe - trcmin)
      coccflag = 0.5 + sign(0.5,biococc - trcmin)
      caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)

!     limit tracers to positive values
      bioin(:) = max(bioin(:), trcmin)
      biopo4 = max(biopo4, trcmin)
      biophyt = max(biophyt, trcmin)
      biozoop = max(biozoop, trcmin)
      biodetr = max(biodetr, trcmin)
      biodic = max(biodic, trcmin)
      biono3 = max(biono3, trcmin)
      biodop = max(biodop, trcmin)
      biodon = max(biodon, trcmin)
      biodiaz = max(biodiaz, trcmin)
      biodin15 = max(biodin15, trcmin)
      biodon15 = max(biodon15, trcmin)
      biophytn15 = max(biophytn15, trcmin)
      biococcn15 = max(biococcn15, trcmin)
      biozoopn15 = max(biozoopn15, trcmin)
      biodetrn15 = max(biodetrn15, trcmin)
      biodiazn15 = max(biodiazn15, trcmin)
      biodic13 = max(biodic13, trcmin)
      biophytc13 = max(biophytc13, trcmin)
      biococcc13 = max(biococcc13, trcmin)
      biocaco3c13 = max(biocaco3c13, trcmin)
      biozoopc13 = max(biozoopc13, trcmin)
      biodetrc13 = max(biodetrc13, trcmin)
      biodoc13 = max(biodoc13, trcmin)
      biodiazc13 = max(biodiazc13, trcmin)
      biococc = max(biococc, trcmin)
      biocaco3 = max(biocaco3, trcmin)
      biodfe = max(biodfe, trcmin)
      biodetrfe = max(biodetrfe, trcmin)

!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      f1 = exp((-kw - kc*(biophyt
     &     + biodiaz
     &     + biococc
     &     )
     &     - kc_c*biocaco3
     &     )*dzt)

      jmax = abio_P*bct*deffe
      gd = jmax*dayfrac
      u1 = max(gl_O/gd,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2

      avej = gd*(phi1 - phi2)/((kw+kc*(biophyt
     &     + biodiaz
     &     + biococc
     &     )
     &     + kc_c*biocaco3
     &       )*dzt)
!     Make the max grazing rate a function of temperature
!     bctz sets an upper limit on the effects of temp on grazing
!     in contrast to phytoplankton growth rates bct, which are unlimited
      gmax = gbio*bctz
      jmax_D = max(0.,abio_P*(bct - 2.6)*deffe_D)*jdiar
!
      gd_D = max(1.e-14,jmax_D*dayfrac)
      u1 = max(gl_D/gd_D,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*(phi1 - phi2)/((kw+kc*(biophyt+biodiaz
     &       + biococc
     &        )
     &       + kc_c*biocaco3
     &         )*dzt)
      jmax_C = abio_C*bct*deffe_C
      gd_C = jmax_C*dayfrac
      u1 = max(gl_C/gd_C,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_C = gd_C*(phi1 - phi2)/((kw+kc*(biophyt
     &       +  biodiaz
     &       +  biococc)
     &       + kc_c*biocaco3
     &       )*dzt)
      nupt = nupt0*bct
      nupt_D = nupt0_D*bct
      nfixout = 0.0
      expoout = 0.0
      grazout = 0.0
      morpout = 0.0
      morzout = 0.0
      graz_Det_out = 0.0
      graz_Zout = 0.0
      rn15expoout = 0.0
      rc13expoout = 0.0
      calproout = 0.0
      nuct = nuct0*bct
      npp_Cout = 0.0
      graz_Cout = 0.0
      morp_Cout = 0.0
      morpt_Cout = 0.0
      calattout = 0.0
      disslout = 0.0
      expocaco3out = 0.0
      expofeout = 0.0
      remifeout = 0.0
      nppout = 0.0
      npp_dopout = 0.0
      morptout = 0.0
      remiout = 0.0
      excrout = 0.0
      npp_Dout = 0.0
      npp_D_dopout = 0.0
      graz_Dout = 0.0
      morpt_Dout = 0.0
      morp_Dout = 0.0
      npp_C_dopout = 0.0

      do n=1,nbio

        p1 = min(biophyt,pmax)
        p2 = max(0.0,biophyt - pmax)
        kfevar = (kfemin*p1+ kfemax*p2)/(p1 + p2)
        deffe = biodfe/(kfevar + biodfe)
        jmax = abio_P*bct*deffe
        p1_C = min(biococc,pmax_C)
        p2_C = max(0.0,biococc - pmax_C)
        kfevar_C = (kfemin_C*p1_C + kfemax_C*p2_C)/(p1_C+p2_C)
        deffe_C = biodfe/(kfevar_C + biodfe)
        jmax_C = abio_C*bct*deffe_C
        deffe_D = biodfe/(kfe_D + biodfe)
        jmax_D = max(0.,abio_P*(bct - 2.6)*deffe_D)*jdiar
!       growth rate of phytoplankton
!       consume DOP when it is more efficient
        limP_dop = hdop*biodop/(k1p_P + biodop)
        limP_po4 = biopo4/(k1p_P + biopo4)
        dopupt_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP = limP_dop*dopupt_flag + limP_po4*(1.-dopupt_flag)
        u_P = min(avej, jmax*limP)
        limP_dop = hdop*biodop/(k1p_C + biodop)
        limP_po4 = biopo4/(k1p_C + biopo4)
        dopupt_C_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP_C = limP_dop*dopupt_C_flag + limP_po4*(1.-dopupt_C_flag)
        u_C = min(avej, jmax_C*limP_C)
        po4P = jmax*biopo4/(k1p_P + biopo4)

!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
        u_C = min(u_C, jmax_C*biono3/(k1n_C + biono3))
        no3P = jmax*biono3/(k1n + biono3)
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        u_D = min(avej_D, jmax_D*limP)
        dopupt_D_flag = dopupt_flag
        po4_D = jmax_D*biopo4/(k1p_P + biopo4)
!       Set the grazing coefficients for the N case
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop
     &         + zprefDiaz*biodiaz + kzoo
     &         + zprefC*biococc
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
        ing_D = zprefDiaz/thetaZ
        ing_C = zprefC/thetaZ
        npp = u_P*biophyt
        npp_C = u_C*biococc
        dopupt = npp*dopupt_flag
        dopupt_C = npp_C*dopupt_C_flag
        npp_D = max(0.,u_D*biodiaz)
!       grazing on diazotrophs
        g_D = gmax*ing_D*biodiaz
        graz_D = g_D*biozoop
        morpt_D = nupt_D*biodiaz ! linear mortality
        morp_D = nup_D*biodiaz*biodiaz
c        no3upt_D = biono3/(k1n + biono3)*npp_D ! nitrate uptake
        no3upt_D = (0.5+0.5*tanh(biono3-5.))*npp_D ! nitrate uptake
        dopupt_D = npp_D*dopupt_D_flag
!       grazing on P
        g_P = gmax*ing_P*biophyt
        graz = g_P*biozoop
!       grazing on Z
        g_Z = gmax*ing_Z*biozoop
        graz_Z = g_Z*biozoop
!       grazing on Detritus
        g_Det = gmax*ing_Det*biodetr
        graz_Det = g_Det*biozoop
!
        morp = nup*biophyt
        morpt = nupt*biophyt
        recy_don = nudon*bct*biodon
        recy_dop = nudop*bct*biodop
        morz = nuz*biozoop*biozoop
        remi = nud*bct*biodetr
        expo = wwd*biodetr
        g_C =  gmax*ing_C*biococc
        graz_C = g_C*biozoop
        morp_C = nuc*biococc
        morpt_C = nuct*biococc
        dissl = biocaco3*dissk1
        expocaco3 = wwc*biocaco3
!   Calculation of average light in mixed layer for calculation of Chl diagnostic
        par = 0.43  ! fraction of photosythetically active radiation
        irrtop=gl/2/par
        kirr=-kw - kc*(bioin(imobiphyt)
     &      + bioin(imobidiaz)
     &      + bioin(imobicocc)
     &       )
     &      - kc_c*bioin(imobicaco3)

        aveirr=-1/dzt/kirr*(irrtop-irrtop*exp(kirr*dzt))
!       remineralization of iron from organic matter
        remife=nud*bct*biodetrfe

!       Scavenging of dissolved iron is based on Honeymoon et al. (1988)
!       and Parekh et al. (2004).
!       o2flag is zero for o2 < o2min and one otherwise
        o2flag = 0.5 + sign(0.5, o2-o2min)
        fepa = (1.0 + kfeleq * (lig - biodfe))*o2flag

        feprime = ((-fepa +(fepa
     &         * fepa + 4.0 * kfeleq
     &         * biodfe)**(0.5)) /(2.0 * kfeleq))*o2flag

        feorgads = (kfeorg*(((biodetr * detrflag
     &         )*mc*redctn)**0.58)*feprime)*o2flag
        fecol = kfecol*feprime*o2flag

        expofe = wwd * biodetrfe

!

        graz = graz*phytflag*zoopflag
     &         *phytn15flag*zoopn15flag
        graz_Z = graz_Z*zoopflag
     &          *zoopn15flag
        graz_Det = graz_Det*detrflag*zoopflag
     &            *detrn15flag*zoopn15flag
        morp = morp*phytflag
     &        *phytn15flag
        morpt = morpt*phytflag
     &        *phytn15flag
        morz = morz*zoopflag
     &        *zoopn15flag
        remi = remi*detrflag
     &        *detrn15flag
        expo = expo*detrflag
     &        *detrn15flag
        recy_dop = recy_dop*dopflag
           npp = npp*no3flag*(dopupt_flag*dopflag
     &         + (1.-dopupt_flag)*po4flag)
     &           *din15flag
           npp_C = npp_C*no3flag*(dopupt_C_flag*dopflag
     &           + (1.-dopupt_C_flag)*po4flag)
     &             *din15flag
           npp_D = npp_D*(dopupt_D_flag*dopflag
     &           + (1.-dopupt_D_flag)*po4flag)
     &             *din15flag
        graz_D = graz_D*diazflag*zoopflag
     &          *diazn15flag*zoopn15flag
        morpt_D = morpt_D*diazflag
     &           *diazn15flag
        morp_D = morp_D*diazflag
     &           *diazn15flag
        no3upt_D = no3upt_D*no3flag
     &            *din15flag
        recy_don = recy_don*donflag
     &           *don15flag
c Andreas commented out (don't think we need these)
c        calpro = calpro*caco3flag*coccflag*zoopflag
        dissl = dissl*caco3flag ! AS add caco3c13flag
        expocaco3 = expocaco3*caco3flag
c        calatt = calatt*caco3flag*coccflag*zoopflag
        graz_C = graz_C*coccflag*zoopflag
     &        *coccn15flag*zoopn15flag
        morp_C = morp_C*coccflag
     &        *coccn15flag
        morpt_C = morpt_C*coccflag
     &        *coccn15flag
        remife=remife*detrfeflag
        feorgads=feorgads*dfeflag
        expofe=expofe*detrfeflag
        fecol=fecol*dfeflag

!     digestion of grazed material
        dig_P = gamma1*graz
        dig_Z = gamma1*graz_Z
        dig_Det = gamma1*graz_Det
        dig_C = gamma1*graz_C
        dig = dig_P + dig_Z + dig_Det
     &      + dig_C
!     excretion based on growth efficiency
        excr_P = gamma1*(1-geZ)*graz
        excr_Z = gamma1*(1-geZ)*graz_Z
        excr_Det = gamma1*(1-geZ)*graz_Det
        excr_C = gamma1*(1-geZ)*graz_C
        excr = excr_P + excr_Z + excr_Det
     &       + excr_C
!     sloppy feeding
        sf_P = (1-gamma1)*graz
        sf_Z = (1-gamma1)*graz_Z
        sf_Det = (1-gamma1)*graz_Det
        sf_C = (1-gamma1)*graz_C
        sf = sf_P + sf_Z + sf_Det
     &       + sf_C
!     digestion of grazed material
        dig_D = gamma1*graz_D*(redntp/diazntp)
        dig = dig + dig_D

!     excretion based on growth efficiency
        excr_D = gamma1*(1-geZ)*graz_D*(redntp/diazntp)
        excr = excr + excr_D
!       excrete "extra" non-Redfield diazotroph N
        nr_excr_D = gamma1*graz_D*(1-(redntp/diazntp))
     &       + (1-gamma1)*graz_D*(1-(redntp/diazntp))

!     sloppy feeding
        sf_D = (1-gamma1)*graz_D*(redntp/diazntp)
        sf = sf + sf_D
!       calculate isotope parameters
!       See Somes et al., 2010, GBC for details/results
        uno3 = npp*dtbio/biono3
        uno3 = min(uno3, 0.999)
        uno3 = max(uno3, trcmin)
        rno3 = biodin15/(biono3-biodin15)
        rno3 = min(rno3, 2*rn15std)
        rno3 = max(rno3, rn15std/2.)
        bassim = rno3 + eps_assim*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
        fcassim = bassim/(1+bassim)

        udon = recy_don*dtbio/biodon
        udon = min(udon, 0.999)
        udon = max(udon, trcmin)
        rdon = biodon15/(biodon-biodon15)
        rdon = min(rdon, 2*rn15std)
        rdon = max(rdon, rn15std/2.)
        brecy = rdon + eps_recy*(1-udon)/udon*log(1-udon)*rdon/1000.
        fcrecy = brecy/(1+brecy)

        rzoop = biozoopn15/(biozoop-biozoopn15)
        rzoop = min(rzoop, 2.*rn15std)
        rzoop = max(rzoop, rn15std/2.)
        bexcr = rzoop - eps_excr*rzoop/1000.
        fcexcr = bexcr/(1+bexcr)

        bnfix = rn15std - eps_nfix*rn15std/1000.
        fcnfix = bnfix/(1+bnfix)

        rtdin15 = biodin15/biono3
        rtdin15 = min(rtdin15, 2*rn15std/(1+rn15std))
        rtdin15 = max(rtdin15, rn15std/(1+rn15std)/2.)

        rtdon15 = biodon15/biodon
        rtdon15 = min(rtdon15, 2*rn15std/(1+rn15std))
        rtdon15 = max(rtdon15, rn15std/(1+rn15std)/2.)

        rtphytn15 = biophytn15/biophyt
        rtphytn15 = min(rtphytn15, 2.*rn15std/(1+rn15std))
        rtphytn15 = max(rtphytn15, rn15std/(1+rn15std)/2.)

        rtcoccn15 = biococcn15/biococc
        rtcoccn15 = min(rtcoccn15, 2.*rn15std/(1+rn15std))
        rtcoccn15 = max(rtcoccn15, rn15std/(1+rn15std)/2.)

        rtzoopn15 = biozoopn15/biozoop
        rtzoopn15 = min(rtzoopn15, 2.*rn15std/(1+rn15std))
        rtzoopn15 = max(rtzoopn15, rn15std/(1+rn15std)/2.)

        rtdetrn15 = biodetrn15/biodetr
        rtdetrn15 = min(rtdetrn15, 2.*rn15std/(1+rn15std))
        rtdetrn15 = max(rtdetrn15, rn15std/(1+rn15std)/2.)

        rtdiazn15 = biodiazn15/biodiaz
        rtdiazn15 = min(rtdiazn15, 2.*rn15std/(1+rn15std))
        rtdiazn15 = max(rtdiazn15, rn15std/(1+rn15std)/2.)

        rdic13 = biodic13/(biodic-biodic13)
        rdic13 = min(rdic13, 2.*rc13std)
        rdic13 = max(rdic13, 0.5*rc13std)
        bc13npp = ac13b*rdic13
        fcnpp = bc13npp/(1+bc13npp)

        rtdic13 = biodic13/biodic
        rtdic13 = min(rtdic13, 2*rc13std/(1+rc13std))
        rtdic13 = max(rtdic13, 0.5*rc13std/(1+rc13std))

        rtphytc13 = biophytc13/(biophyt*redctn)
        rtphytc13 = min(rtphytc13, 2.*rc13std/(1+rc13std))
        rtphytc13 = max(rtphytc13, 0.5*rc13std/(1+rc13std))

        rtcoccc13 = biococcc13/(biococc*redctn)
        rtcoccc13 = min(rtcoccc13, 2.*rc13std/(1+rc13std))
        rtcoccc13 = max(rtcoccc13, 0.5*rc13std/(1+rc13std))

        rtcaco3c13 = biocaco3c13/biocaco3
        rtcaco3c13 = min(rtcaco3c13, 2.*rc13std/(1+rc13std))
        rtcaco3c13 = max(rtcaco3c13, 0.5*rc13std/(1+rc13std))

        rtzoopc13 = biozoopc13/(biozoop*redctn)
        rtzoopc13 = min(rtzoopc13, 2.*rc13std/(1+rc13std))
        rtzoopc13 = max(rtzoopc13, 0.5*rc13std/(1+rc13std))

        rtdetrc13 = biodetrc13/(biodetr*redctn)
        rtdetrc13 = min(rtdetrc13, 2.*rc13std/(1+rc13std))
        rtdetrc13 = max(rtdetrc13, 0.5*rc13std/(1+rc13std))

        rtdoc13 = biodoc13/(biodon*redctn)
        rtdoc13 = min(rtdoc13, 2*rc13std/(1+rc13std))
        rtdoc13 = max(rtdoc13, 0.5*rc13std/(1+rc13std))

        rtdiazc13 = biodiazc13/(biodiaz*redctn)
        rtdiazc13 = min(rtdiazc13, 2.*rc13std/(1+rc13std))
        rtdiazc13 = max(rtdiazc13, 0.5*rc13std/(1+rc13std))

!  net formation of attached (living) tests
!     total primary production by coccs
        calatt = (npp_C - morpt_C - morp_C - graz_C
! total growth of zooplankton
     &         + dig - morz - graz_Z - excr
! convert to carbon units mmol C
     &          )*capr*redctn*1.e3

!   formation of detached (dead) tests, or PIC
        calpro = (sf_C + sf_Z + morp_C + morz)*capr*redctn*1.e3 !stay in mmol
!       nutrients equation
        biopo4 = biopo4 + dtbio*(redptn*(excr + remi
     &           + (1.-dfrt)*morpt - (npp-dopupt)
     &           + (1.-dfrt)*morpt_C - (npp_C-dopupt_C)
     &            )
     &         + diazptn*(morpt_D - (npp_D-dopupt_D)) + recy_dop
     &          )

!       DOP equation
        biodop = biodop + dtbio*(redptn*(
     &             dfr*morp + dfrt*morpt - dopupt
     &           + dfr*morp_C + dfrt*morpt_C - dopupt_C
     &            )
     &         - diazptn*dopupt_D - recy_dop)

!       phytoplankton equation
        biophyt = biophyt + dtbio*(npp - morp - graz - morpt)
!       zooplankton equation
        biozoop = biozoop + dtbio*(dig - morz - graz_Z - excr)

!       detritus equation
        biodetr = biodetr + dtbio*((1.-dfr)*morp + sf + morz - remi
     &          - graz_Det - expo + impo + morp_D*(redntp/diazntp)
     &          + (1.-dfr)*morp_C
     &           )

!       DIC equation
        biodic = biodic + dtbio*redctn*(excr + remi
     &         + (1.-dfrt)*morpt - npp
     &         + (1.-dfrt)*morpt_C - npp_C
     &         + morpt_D - npp_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )
!       nitrate (NO3) equation
        biono3 = biono3 + dtbio*(excr + remi
     &         + (1.-dfrt)*morpt - npp
     &         + (1.-dfrt)*morpt_C - npp_C
     &         + morpt_D - no3upt_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )
!       DON equation
        biodon = biodon + dtbio*(dfr*morp + dfrt*morpt - recy_don
     &         + dfr*morp_C + dfrt*morpt_C
     &       )
!       diazotroph equation
        biodiaz = biodiaz + dtbio*(npp_D - morp_D - morpt_D - graz_D)
!       coccolithophores equation
        biococc = biococc + dtbio*(npp_C - morp_C - graz_C
     &       - morpt_C)
!     calcite equation
c        biocaco3 = biocaco3 - dissl + dtbio*(calpro-expocaco3+impocaco3)
c     Andreas put dissl in bracket
        biocaco3 = biocaco3 + dtbio*(calpro-dissl-expocaco3+impocaco3)
!       dissolved iron equation
        biodfe = biodfe + dtbio*(rfeton*(excr + (1.-dfrt)*morpt
     &            - npp + morpt_D - npp_D + recy_don + nr_excr_D
     &            + morp_D*(1-(redntp/diazntp))) - feorgads
     &            + remife - fecol
     &            + rfeton*((1.-dfrt)*morpt_C - npp_C)
     &            )
!     particulate iron equation
        biodetrfe = biodetrfe + dtbio*(rfeton*(sf + (1.-dfr)*morp
     &            + morp_D*(redntp/diazntp) + morz - graz_Det)
     &            + feorgads + fecol - remife - expofe + impofe
     &            + rfeton*(1.-dfr)*morp_C
     &       )
!       isotope equations
        biodin15 = biodin15 + dtbio*(rtphytn15*(1.-dfrt)*morpt
     &       + rtcoccn15*(1.-dfrt)*morpt_C - fcassim*npp_C
     &       + fcexcr*excr + rtdiazn15*morpt_D + rtdiazn15*nr_excr_D
     &       + rtdiazn15*morp_D*(1-(redntp/diazntp)) + rtdetrn15*remi
     &       + fcrecy*recy_don - fcassim*npp - fcassim*no3upt_D)

        biodon15 = biodon15 + dtbio*(dfr*rtphytn15*morp
     &       + dfr*rtcoccn15*morp_C + dfrt*rtcoccn15*morpt_C
     &       + dfrt*rtphytn15*morpt - fcrecy*recy_don)

        biophytn15 = biophytn15 + dtbio*(fcassim*npp - rtphytn15*morp
     &       - rtphytn15*graz - rtphytn15*morpt)
        biococcn15 = biococcn15 + dtbio*(fcassim*npp_C-rtcoccn15*morp_C
     &       - rtcoccn15*graz_C - rtcoccn15*morpt_C)
        biozoopn15 = biozoopn15 + dtbio*(rtphytn15*dig_P
     &       + rtcoccn15*dig_C
     &       + rtzoopn15*dig_Z + rtdetrn15*dig_Det + rtdiazn15*dig_D
     &       - rtzoopn15*morz - rtzoopn15*graz_Z - fcexcr*excr)

        biodetrn15 = biodetrn15 + dtbio*(rtphytn15*(1.-dfr)*morp
     &       + rtcoccn15*(1.-dfr)*morp_C + rtcoccn15*sf_C
     &       + rtphytn15*sf_P + rtzoopn15*sf_Z + rtdetrn15*sf_Det
     &       + rtdiazn15*sf_D + rtzoopn15*morz - rtdetrn15*remi
     &       - rtdetrn15*graz_Det - rtdetrn15*expo + rn15impo*impo
     &       + rtdiazn15*morp_D*(redntp/diazntp))

        biodiazn15 = biodiazn15 + dtbio*(fcnfix*(npp_D-no3upt_D)
     &       + fcassim*no3upt_D - rtdiazn15*morp_D - rtdiazn15*graz_D
     &       - rtdiazn15*morpt_D)
!     !!!!!!!!!!!!!!!!!! isotope equations
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*(1.-dfrt)*morpt
     &       + rtzoopc13*excr + rtdiazc13*morpt_D + rtdiazc13*nr_excr_D
     &       + rtdiazc13*morp_D*(1-(redntp/diazntp)) + rtdetrc13*remi
     &       + rtcoccc13*(1.-dfrt)*morpt_C - fcnpp*npp_C
     &       + rtdoc13*recy_don - fcnpp*npp - fcnpp*npp_D)

        biodoc13 = biodoc13 + dtbio*redctn*(dfr*rtphytc13*morp
     &       + rtcoccc13*(dfr*morp_C + dfrt*morpt_C)
     &       + rtphytc13*dfrt*morpt - rtdoc13*recy_don)

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp
     &       - rtphytc13*morp - rtphytc13*graz - rtphytc13*morpt)

        biococcc13 = biococcc13 + dtbio*redctn*(fcnpp*npp_C
     &       - rtcoccc13*morp_C - rtcoccc13*graz_C - rtcoccc13*morpt_C)
        biocaco3c13 = biocaco3c13 + dtbio*(rtdic13*calpro
     &       - rtcaco3c13*dissl - rtcaco3c13*expocaco3
     &                + rcaco3c13impo*impocaco3)

        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P
     &      + rtcoccc13*dig_C
     &      + rtzoopc13*dig_Z + rtdetrc13*dig_Det + rtdiazc13*dig_D
     &       - rtzoopc13*morz - rtzoopc13*graz_Z - rtzoopc13*excr)

        biodetrc13 = biodetrc13 + dtbio*redctn*(rtphytc13*(1.-dfr)*morp
     &       + rtcoccc13*(1.-dfr)*morp_C + rtcoccc13*sf_C
     &       + rtphytc13*sf_P + rtzoopc13*sf_Z + rtdetrc13*sf_Det
     &       + rtzoopc13*morz - rtdetrc13*remi
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo*impo
     &       + rtdiazc13*morp_D*(redntp/diazntp))

        biodiazc13 = biodiazc13 + dtbio*redctn*(fcnpp*npp_D
     &       - rtdiazc13*(morp_D + graz_D + morpt_D))
        expoout = expoout + expo
        rn15expoout = rn15expoout + rtdetrn15
        rc13expoout = rc13expoout + rtdetrc13
        rcaco3c13expoout = rcaco3c13expoout + rtcaco3c13
        calproout = calproout + calpro
        nfixout = nfixout + npp_D - no3upt_D
        expofeout = expofeout + expofe
        remifeout = remifeout + remife
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
        graz_Det_out = graz_Det_out + graz_Det
        graz_Zout = graz_Zout + graz_Z
        nppout = nppout + npp
        morptout = morptout + morpt
        remiout = remiout + remi
        excrout = excrout + excr
        npp_dopout = npp_dopout + npp*dopupt_flag
        npp_C_dopout = npp_C_dopout + npp_C*dopupt_C_flag
        npp_Dout = npp_Dout + npp_D
        npp_D_dopout = npp_D_dopout + npp_D*dopupt_D_flag
        graz_Dout = graz_Dout + graz_D
        morpt_Dout = morpt_Dout + morpt_D
        morp_Dout = morp_Dout + morp_D
        calattout = calattout + calatt
        disslout = disslout + dissl
        expocaco3out = expocaco3out + expocaco3
        npp_Cout = npp_Cout + npp_C
        graz_Cout = graz_Cout + graz_C
        morp_Cout = morp_Cout + morp_C
        morpt_Cout = morpt_Cout + morpt_C
!     flags prevent negative values by setting outgoing fluxes to
!     zero if tracers are lower than trcmin
        if (po4flag .eq. 1) po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
        if (phytflag .eq. 1) phytflag = 0.5 + sign(0.5,biophyt - trcmin)
        if (zoopflag .eq. 1) zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
        if (detrflag .eq. 1) detrflag = 0.5 + sign(0.5,biodetr - trcmin)
        if (no3flag .eq. 1) no3flag = 0.5 + sign(0.5,biono3 - trcmin)
        if (dopflag .eq. 1) dopflag = 0.5 + sign(0.5,biodop - trcmin)
        if (donflag .eq. 1) donflag = 0.5 + sign(0.5,biodon - trcmin)
        if (diazflag .eq. 1) diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
        if (din15flag .eq. 1)
     &     din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
        if (don15flag .eq. 1)
     &     don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
        if (phytn15flag .eq. 1)
     &     phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
        if (coccn15flag .eq. 1)
     &     coccn15flag = 0.5 + sign(0.5, biococcn15 - trcmin)
        if (zoopn15flag .eq. 1)
     &     zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
        if (detrn15flag .eq. 1)
     &     detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
        if (diazn15flag .eq. 1)
     &     diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
        if (caco3flag .eq. 1)
     &     caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)
        if (coccflag .eq. 1)
     &     coccflag = 0.5 + sign(0.5,biococc - trcmin)
        if (dfeflag .eq. 1) dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
        if (detrfeflag .eq. 1) detrfeflag = 0.5
     &                                 + sign(0.5,biodetrfe - trcmin)
        if (dic13flag .eq. 1)
     &       dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
        if (phytc13flag .eq. 1)
     &       phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
        if (coccc13flag .eq. 1)
     &       coccc13flag = 0.5 + sign(0.5,biococcc13 - trcmin)
        if (zoopc13flag .eq. 1)
     &       zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
        if (detrc13flag .eq. 1)
     &       detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
        if (doc13flag .eq. 1)
     &       doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
        if (diazc13flag .eq. 1)
     &       diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
      enddo

      bioout(imobipo4) = biopo4 - bioin(imobipo4)
      bioout(imobiphyt) = biophyt - bioin(imobiphyt)
      bioout(imobizoop) = biozoop - bioin(imobizoop)
      bioout(imobidetr) = biodetr - bioin(imobidetr)
      bioout(imobidic) = biodic - bioin(imobidic)
      bioout(imobidop) = biodop - bioin(imobidop)
      bioout(imobino3) = biono3 - bioin(imobino3)
      bioout(imobidon) = biodon - bioin(imobidon)
      bioout(imobidiaz) = biodiaz - bioin(imobidiaz)
      bioout(imobidin15) = biodin15 - bioin(imobidin15)
      bioout(imobidon15) = biodon15 - bioin(imobidon15)
      bioout(imobiphytn15) = biophytn15 - bioin(imobiphytn15)
      bioout(imobizoopn15) = biozoopn15 - bioin(imobizoopn15)
      bioout(imobidetrn15) = biodetrn15 - bioin(imobidetrn15)
      bioout(imobidiazn15) = biodiazn15 - bioin(imobidiazn15)
      bioout(imobicoccn15) = biococcn15 - bioin(imobicoccn15)
      bioout(imobicocc) = biococc - bioin(imobicocc)
      bioout(imobicaco3) = biocaco3 - bioin(imobicaco3)
      bioout(imobidfe) = biodfe - bioin(imobidfe)
      bioout(imobidetrfe) = biodetrfe - bioin(imobidetrfe)
      bioout(imobidic13) = biodic13 - bioin(imobidic13)
      bioout(imobiphytc13) = biophytc13 - bioin(imobiphytc13)
      bioout(imobizoopc13) = biozoopc13 - bioin(imobizoopc13)
      bioout(imobidetrc13) = biodetrc13 - bioin(imobidetrc13)
      bioout(imobidoc13) = biodoc13 - bioin(imobidoc13)
      bioout(imobidiazc13) = biodiazc13 - bioin(imobidiazc13)
      bioout(imobicoccc13) = biococcc13 - bioin(imobicoccc13)
      bioout(imobicaco3c13) = biocaco3c13 - bioin(imobicaco3c13)
      return
      end
