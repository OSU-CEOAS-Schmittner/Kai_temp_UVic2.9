! source file: /raid24/aschmitt/UVic2.9/MOBI1.9/nobio/updates/npzd.h
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
     &                     )
      common /npzd_i/ nbio

      integer imobipo4, imobiphyt, imobizoop, imobidetr
      common /npzd_i/ imobipo4, imobiphyt, imobizoop, imobidetr

      real trcmin
      parameter (trcmin=5e-12)

