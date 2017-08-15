! source file: /data/home/kai/dev/UVic2.9/updates/npzd.h
!====================== include file "npzd.h" =========================

!   variables for npzd model

!   ntnpzd   = number of npzd tracers
!   nbio     = number of npzd timesteps per ocean timestep
!   trcmin   = minimum tracer for flux calculations
!   tap      = 2*alpha*par with
!   alpha    = initial slope P-I curve [(W/m^2)^(-1)/day] and
!   par      = fraction of photosythetically active radiation
!   kw       = light attenuation due to water [1/m]
!   kc       = light attenuation by phytoplankton [1/(m*mmol m-3)]
!   ki       = light attenuation through sea ice & snow
!   abio     = maximum growth rate parameter [1/day]
!   bbio     = b
!   cbio     = [1/deg_C]
!   k1n      = half saturation constant for N uptake [mmol m-3]
!   nup      = specific mortality rate (Phytoplankton) [day-1]
!   gamma1   = assimilation efficiency (zpk)
!   gbio     = maximum grazing rate at 0 deg C [day-1]
!   nuz      = quadratic mortality (zpk)
!   nud0     = remineralization rate [day-1]
!   nudop0   = DON remineralization rate [day-1]
!   nudon0   = DOP remineralization rate [day-1]
!   wd       = sinking speed of detritus [m day-1]
!   ztt      = depth to top of grid cell [cm]
!   rkwz     = reciprical of light attenuation times grid depth
!   dtnpzd   = time step of biology
!   capr     = carbonate to carbon production ratio
!   dcaco3   = remineralisation depth of calcite [cm]
!   rcak     = array used in calculating calcite remineralization
!   rcab     = array used in calculating bottom calcite remineralization
!   nupt0    = specific mortality rate (Phytoplankton) [1/day]
!   wd0      = sinking speed of detritus at surface [m/day]
!   mw       = sinking speed increase with depth [1/day]
!   mw_c       = calcite sinking speed increase with depth [1/day]
!   mwz      = sinking speed increase depth cut-off (cm)
!   k1p_P    = half saturation constant for P uptake phytoplankton
!   jdiar    = factor reducing the growth rate of diazotrophs
!   redctn   = C/N Redfield ratio (includes mol to mmol conversion)
!   redctp   = C/P Redfield ratio (includes mol to mmol conversion)
!   redptn   = P/N Redfield ratio
!   redntp   = N/P Redfield ratio
!   redotn   = O/N Redfield ratio (includes mol to mmol conversion)
!   redotp   = O/P Redfield ratio (includes mol to mmol conversion)
!   rnbio    = reciprical of nbio
!   rdtts    = reciprical of dtts [s-1]
!   dtbio    = npzd time step [s]
!   rnpp     = rate of net primary production [nmol cm-3 s-1]
!   rgraz    = rate of grazing [nmol cm-3 s-1]
!   rmorp    = rate of mortality of phytoplankton [nmol cm-3 s-1]
!   rmorz    = rate of mortality of zooplankton [nmol cm-3 s-1]
!   rremi    = rate of remineralization [nmol cm-3 s-1]
!   rexcr    = rate of excretion [nmol cm-3 s-1]
!   rexpo    = rate of export through the bottom [nmol cm-3 s-1]
!   rnpp_D   = npp for diazotraphs [nmol cm-3 s-1]
!   rgraz_D  = rgraz for diazotraphs [nmol cm-3 s-1]
!   rmorpt_D = rmorp for diazotraphs [nmol cm-3 s-1]
!   rnfix    = rate of nitrogen fixation [nmol cm-3 s-1]
!   rdeni    = rate of denitrification [nmol cm-3 s-1]
!   kzoo     = half saturation constant for Z grazing
!   zprefP   = Z preference for grazing on P
!   zprefDiaz   = Z preference for grazing on Diaz
!   zprefZ   = Z preference for grazing on other Z
!   zprefDet = Z preference for grazing on Detritus
!   rgraz_Det = rate of grazing on Detritus [nmol cm-3 s-1]
!   rgraz_Z   = rate of grazing on other Zooplankton [nmol cm-3 s-1]
!   geZ      = growth efficiency of zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! New diagnostic output
!   ravej      = light-dependant growth rate of P
!   ravej_D    = light-dependant growth rate of Diaz
!   rgmax      = temp-dependant growth rate of zoo
!   rno3P      = nitrate-dependant growth rate of P
!   rpo4P       = phosphate-dependant growth rate of P
!   rpo4_D     = phosphate-dependant growth rate of D
!
!   fe_dissolved = dissolved iron concentration
!   kfe = Fe limitation half saturation parameter
!   kfe_D = Fe limitation half sat. param. for diaz.
!
!++++ Climate engineering ++++++++++++++
!   kpipe = ocean pipe coordinate
!   kpipe_fe = ocean pipe coordinate for fe limitation
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Dynamic iron cycle
!   kfeleq = Fe-ligand stability constant [m^3/(mmol ligand)]
!   lig = Ligand concentration  [mmol/m^3]
!   thetamaxhi = Maximum Chl:C ratio, abundant iron [gChl/(gC)]
!   thetamaxlo = Maximum Chl:C ratio, extreme iron limitation [gChl/(gC)]
!   alphamax = Maximum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   alphamin = Minimum intital slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   mc = Molar mass of carbon [g/mol]
!   fetopsed = Fe:P ratio for sedimentary iron source [molFe/molP]
!   o2min = Minimum O2 concentration for aerobic respiration [mmolO_2/m^3]
!   kfeorg = Organic-matter dependent scavenging rate [(m^3/(gC s))^0.58]
!   kfecol = Colloidal production and precipitation rate [s^-1]

      integer ntnpzd, nbio
      parameter (ntnpzd = 4 ! po4, phyt, zoop, detr
     &                  + 1 ! dic
     &                  + 4 ! dic13, phytc13, zoopc13, detrc13
     &                  +2  ! doc13, diazc13
     &                  +2  ! coccc13, caco3c13
     &                  + 4 ! no3, diaz, don, dop
     &                  + 6 ! din15, don15, phytn15, zoopn15, detrn15, diazn15
     &                  + 1 ! coccn15
     &                   +2 ! cocc, caco3
     &                  + 2 ! dfe, detrfe
     &                     )
      common /npzd_i/ nbio

      integer imobipo4, imobiphyt, imobizoop, imobidetr
      common /npzd_i/ imobipo4, imobiphyt, imobizoop, imobidetr
      integer imobidic
      common /npzd_i/ imobidic
      integer imobidic13, imobidoc13, imobiphytc13
      integer imobizoopc13, imobidetrc13, imobidiazc13
      common /npzd_i/ imobidic13, imobidoc13, imobiphytc13
      common /npzd_i/ imobizoopc13, imobidetrc13, imobidiazc13
      integer imobicoccc13
      common /npzd_i/ imobicoccc13
      integer imobicaco3c13
      common /npzd_i/ imobicaco3c13
      integer imobidop, imobino3, imobidon, imobidiaz
      common /npzd_i/ imobidop, imobino3, imobidon, imobidiaz
      integer imobidin15, imobidon15, imobiphytn15
      integer imobizoopn15, imobidetrn15, imobidiazn15
      common /npzd_i/ imobidin15, imobidon15, imobiphytn15
      common /npzd_i/ imobizoopn15, imobidetrn15, imobidiazn15
      integer imobicoccn15
      common /npzd_i/ imobicoccn15
      integer imobicocc
      common /npzd_i/ imobicocc
      integer imobicaco3
      common /npzd_i/ imobicaco3
      integer imobidfe, imobidetrfe
      common /npzd_i/ imobidfe, imobidetrfe

      real trcmin
      parameter (trcmin=5e-12)

      real rn15std

      parameter (rn15std=0.0036765)
      real rc13std
      parameter (rc13std=0.0112372)

      real tap, kw, kc, ki, abio_P, bbio, cbio, k1n, nup, gamma1, gbio
      real epsbio, nuz, nud0, LFe, wd, ztt, rkwz, dtnpzd
      real capr, dcaco3, rcak, rcab, nupt0, wd0, k1p_P, jdiar, redctn
      real redctp, redptn, redntp, redotn, redotp, rnbio, rdtts, dtbio
      real rnpp, rgraz, rmorp, rmorpt, rmorz, rremi, rexcr, rexpo
      real rnpp_D, rgraz_D, rmorpt_D, rnfix, kzoo, zprefP, rmorp_D
      real zprefDiaz, zprefZ, zprefDet, rgraz_Det, rgraz_Z, geZ, kfe
      real ravej, ravej_D, rgmax, rno3P, rpo4P, rpo4_D, kfe_D, kpipe
      real kpipe_fe, rwcdeni, rbdeni, rsedrr, sgbdfac, nupt0_D
      real diazntp, diazptn, nup_D, dfr, redotc, redntc, dfrt
      real redptc, rprca, nudop0, nudon0, eps_recy, hdop
      real eps_assim, eps_excr, eps_nfix, eps_wcdeni, eps_bdeni0
      real mw, mwz, mw_c, rnpp_dop, rnpp_D_dop, rnpp_C_dop
      real kfemax, kfemin, pmax
      real kfemax_C, kfemin_C, pmax_C
      real rgraz_Det_B, rexpo_B,rremi_B,bapr

      real abio_C , k1n_C, k1p_C, rnpp_C, nuc, nuct0
      real rgraz_C, rmorp_C, rmorpt_C, tap_C
      real zprefC, kfe_C

      real rcalpro, kcal, wc, wc0, dissk0, rdissl
      real rexpocaco3, rcalatt, rimpocaco3
      real kc_c
!			 , romca, rco3, rco3_sat, rdel_sat

      common /npzd_r/ bapr

      common /npzd_r/ wc0, kcal, dissk0, wc(km)
      common /npzd_r/ kc_c
      common /npzd_r/ rcalpro(kpzd)
!      common /npzd_r/ rcaldiss(kpzd)
      common /npzd_r/ rcalatt(kpzd)
      common /npzd_r/ rdissl(km)
      common /npzd_r/ rexpocaco3(km)
      common /npzd_r/ rimpocaco3(km)
!      common /npzd_r/ romca(km)
!      common /npzd_r/ rco3(km)
!      common /npzd_r/ rco3_sat(km)
!      common /npzd_r/ rdel_sat(km)

      real kfeleq, lig, thetamaxhi, alphamax, alphamin
      real thetamaxlo, mc, fetopsed, o2min, kfeorg, rfeton
      real kfecol, rdeffe_C
      real rfeorgads, rdeffe, rremife, rexpofe, rfeprime
      real rfesed, rbfe, rfecol

      common /npzd_r/ abio_C, k1n_C, k1p_C, nuc, nuct0
      common /npzd_r/ kfe_C, tap_C, mw_c, zprefC
      common /npzd_r/ pmax_C, kfemin_C, kfemax_C
      common /npzd_r/ rnpp_C(kpzd)
      common /npzd_r/ rgraz_C(kpzd)
      common /npzd_r/ rmorp_C(kpzd)
      common /npzd_r/ rmorpt_C(kpzd)

      common /npzd_r/ tap, kw, kc, ki, abio_P, bbio, cbio, k1n, nup
      common /npzd_r/ gamma1, gbio, epsbio, nuz, nud0, LFe, dfr
      common /npzd_r/ wd(km), ztt(km), rkwz(km), dtnpzd, capr
      common /npzd_r/ dcaco3, rcak(km), rcab(km), nupt0, wd0, k1p_P
      common /npzd_r/ jdiar, redctn, redctp, redptn, redntp, redotn
      common /npzd_r/ redotp, rnbio, rdtts, dtbio, geZ
      common /npzd_r/ kzoo, zprefP, zprefDiaz, zprefZ, zprefDet
      common /npzd_r/ kfe, kfe_D
      common /npzd_r/ sgbdfac, nupt0_D, diazntp, diazptn, nup_D
      common /npzd_r/ redotc, redntc, dfrt, redptc, nudop0, nudon0
      common /npzd_r/ eps_assim, eps_excr, eps_nfix, eps_wcdeni
      common /npzd_r/ eps_bdeni0, eps_recy, hdop, mw, mwz
      common /npzd_r/ kfeleq, alphamax, alphamin
      common /npzd_r/ thetamaxhi, thetamaxlo, lig, fetopsed, o2min
      common /npzd_r/ mc, kfeorg, rfeton, kfecol
      common /npzd_r/ kfemax, kfemin, pmax
      common /npzd_r/ rremife(km)
      real fe_hydr
      common /fe_hydr/ fe_hydr(imt,jmt,km)
      common /npzd_r/ rnpp(kpzd), rgraz(kpzd), rmorp(kpzd), rmorpt(kpzd)
      common /npzd_r/ rmorz(kpzd), rexcr(kpzd), rremi(km), rexpo(km)
      common /npzd_r/ rgraz_Det(kpzd), rgraz_Z(kpzd), rsedrr, rprca
      common /npzd_r/ rnpp_dop(kpzd), rnpp_C_dop(kpzd)
      common /npzd_r/ rexpofe(km)
      common /npzd_r/ rnpp_D(kpzd), rgraz_D(kpzd), rmorp_D(kpzd)
      common /npzd_r/ rmorpt_D(kpzd), rnfix(kpzd), rwcdeni(km)
      common /npzd_r/ rbdeni(km), rnpp_D_dop(kpzd)

