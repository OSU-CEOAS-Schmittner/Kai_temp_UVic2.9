! source file: /raid24/aho/UVic2.9/default_comb2/nobio/updates/embm_tsi.F
      subroutine embm_tsi_def (fname, calendar, expnam, runstamp)

!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou

      real c0, c1, c100, c500, c1e3, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

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
      id(1) = id_time

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'

     &, 'time', 'time', 'years since 0-1-1')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_avgper', iou, 1, id, c0, c0, ' ', 'F'
     &,   'averaging period', ' ','day')
      call defvar ('A_sat', iou, 1, id, -c100, c500, ' ', 'F'
     &,   'global average surface air temperature', ' ', 'C')
      call defvar ('A_shum', iou, 1, id, -c100, c100, ' ', 'F'
     &,   'global average surface specific humidity', ' ', '1')
      call defvar ('F_precip', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_evap', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average evaporation', ' ','kg m-2 s-1')
      call defvar ('A_co2', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global average CO2 concentration', ' ','ppm')
      call defvar ('F_co2emit', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global total CO2 emissions', ' ','kg s-1')
      call defvar ('A_maxiter', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'maximum atmosphere solver iterations', ' ','1')
      call defvar ('O_snovol', iou, 1, id, c0, c1e20, ' '
     &,   'F', 'global snow volume', ' ', 'm3')
      call defvar ('O_icevol', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice volume', ' ', 'm3')
      call defvar ('O_icearea', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice area', ' ', 'm2')
      call defvar ('A_satN', iou, 1, id, -c100, c500, ' '
     &,   'F', 'NH average surface air temperature', ' ', 'C')
      call defvar ('A_satS', iou, 1, id, -c100, c500, ' '
     &,   'F', 'SH average surface air temperature', ' ', 'C')
      call defvar ('A_shumN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average surface specific humidity', ' ', '1')
      call defvar ('A_shumS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average surface specific humidity', ' ', '1')
      call defvar ('F_precipN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_precipS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_evapN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average evaporation', ' ','kg m-2 s-1')
      call defvar ('F_evapS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average evaporation', ' ','kg m-2 s-1')
      call defvar ('O_snovolN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH snow volume', ' ', 'm3')
      call defvar ('O_snovolS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH snow volume', ' ', 'm3')
      call defvar ('O_icevolN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice volume', ' ', 'm3')
      call defvar ('O_icevolS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice volume', ' ', 'm3')
      call defvar ('O_iceareaN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice area', ' ', 'm2')
      call defvar ('O_iceareaS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice area', ' ', 'm2')
      call defvar ('A_satL', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over land', ' ', 'C')
      call defvar ('A_satO', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over ocean', ' ', 'C')
      call defvar ('F_precipL', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average precipitation over land', ' ','kg m-2 s-1')
      call defvar ('F_precipO', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average precipitation over ocean', ' ','kg m-2 s-1')
      call defvar ('F_evapL', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average evaporation over land', ' ','kg m-2 s-1')
      call defvar ('F_evapO', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average evaporation over ocean', ' ','kg m-2 s-1')

      call defvar ('F_solins', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'incoming solar insolation', ' ', 'W m-2')
      call defvar ('F_upsens', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'surface upward sensible heat', ' ', 'W m-2')
      call defvar ('F_uplwr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'surface net upward longwave', ' ', 'W m-2')
      call defvar ('F_outlwr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'TOA outgoing longwave', ' ', 'W m-2')
      call defvar ('F_dnswr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net surface downward shortwave (abs.)', ' ', 'W m-2')
      call defvar ('F_absswr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net absorbed shortwave radiation', ' ', 'W m-2')
      call defvar ('F_netrad', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net top of the atmosphere radiation', ' ', 'W m-2')
      call defvar ('A_albplt', iou, 1, id, c0, c1, ' '
     &,    'F', 'plantary albedo', ' ', '1')
      call defvar ('A_albatm', iou, 1, id, c0, c1, ' '
     &,    'F', 'atmospheric albedo', ' ', '1')
      call defvar ('A_albsur', iou, 1, id, c0, c1, ' '
     &,    'F', 'surface albedo', ' ', '1')
      call defvar ('A_albsurL', iou, 1, id, c0, c1, ' '
     &,    'F', 'land surface albedo', ' ', '1')
      call defvar ('A_albsurO', iou, 1, id, c0, c1, ' '
     &,    'F', 'ocean surface albedo', ' ', '1')
      call defvar ('O_tempsur', iou, 1, id, -c100, c500, ' '
     &,   'F', 'global average sea surface temperature', ' ','C')
      call defvar ('O_salsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface salinity', ' ','psu')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine embm_tsi_out (fname, avgper, time, stamp, sat, shum
     &,                        precip, evap, v_oice, a_oice, v_snow
     &,                        v_lice, a_lice, co2ccn, co2emit, dc14ccn
     &,                        dc13ccn
     &,                        cfc11ccn, cfc12ccn, scan, nsat, ssat
     &,                        nshum, sshum, nprecip, sprecip, nevap
     &,                        sevap, v_noice, v_soice, a_noice
     &,                        a_soice, v_nsnow, v_ssnow, v_nlice
     &,                        v_slice, a_nlice, a_slice, lsat, osat
     &,                        lprecip, oprecip, levap, oevap, solins
     &,                        upsens, uplwr, outlwr, dnswr, absswr
     &,                        netrad, palb, aalb, salb, lsalb, osalb
     &,                        sst, sss, ssdic, ssdic13, ssc14, ssalk
     &,                        sso2, sspo4, ssdop, ssno3, ssdon, ssdin15
     &,                        ssdon15, ssdoc13, sscfc11, sscfc12, sulph
     &,                        volc, agg, catm, carbemit, ntrec, ssdfe)
!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname    = file name
!     avgper   = length of averaging period
!     time     = time in years
!     stamp    = time stamp
!     sat, ... = data to be written

!   outputs:
!     ntrec    = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real sat, shum, precip, evap, v_oice, a_oice, v_snow, v_lice
      real a_lice, co2ccn, co2emit, scan, dc14ccn, cfc11ccn, cfc12ccn
      real dc13ccn
      real nsat, ssat, nshum, sshum, nprecip, sprecip, nevap, sevap
      real v_noice, v_soice, a_noice, a_soice, v_nsnow, v_ssnow
      real v_nlice, v_slice, a_nlice, a_slice, lsat, osat, lprecip
      real oprecip, levap, oevap, avgper, solins, upsens, uplwr
      real outlwr, dnswr, absswr, netrad, palb, aalb, salb, lsalb
      real osalb, sst, sss, ssdic, ssc14, ssalk, sso2, sspo4, ssno3
      real sscfc11, sscfc12, sulph, volc, agg, catm, carbemit, time, tmp
      real c0, c1, c100, c1e3, c1e4, c1e6, p1, p001, p035, cal2J, C2K
      real kgsPgyr, ssdop, ssdon, ssdin15, ssdon15, ssdic13, ssdoc13
      real ssdfe

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      C2K = 273.15
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      cal2J = 2.389e-05
      kgsPgyr = 1.e12/(86400.*365.)

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, time, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
      call putvars ('A_sat', iou, ntrec, sat, c1, c0)
      call putvars ('A_shum', iou, ntrec, shum, c1, c0)
      call putvars ('F_precip', iou, ntrec, precip, p1, c0)
      call putvars ('F_evap', iou, ntrec, evap, p1, c0)
      call putvars ('A_co2', iou, ntrec, co2ccn, c1, c0)
      call putvars ('F_co2emit', iou, ntrec, co2emit, c1, c0)
      call putvars ('A_maxiter', iou, ntrec, scan, c1, c0)
      call putvars ('O_icevol', iou, ntrec, v_oice, c1e6, c0)
      call putvars ('O_icearea', iou, ntrec, a_oice, c1e4, c0)
      call putvars ('O_snovol', iou, ntrec, v_snow, c1e6, c0)
      call putvars ('A_satN', iou, ntrec, nsat, c1, c0)
      call putvars ('A_satS', iou, ntrec, ssat, c1, c0)
      call putvars ('A_shumN', iou, ntrec, nshum, c1, c0)
      call putvars ('A_shumS', iou, ntrec, sshum, c1, c0)
      call putvars ('F_precipN', iou, ntrec, nprecip, p1, c0)
      call putvars ('F_precipS', iou, ntrec, sprecip, p1, c0)
      call putvars ('F_evapN', iou, ntrec, nevap, p1, c0)
      call putvars ('F_evapS', iou, ntrec, sevap, p1, c0)
      call putvars ('O_icevolN', iou, ntrec, v_noice, c1e6, c0)
      call putvars ('O_icevolS', iou, ntrec, v_soice, c1e6, c0)
      call putvars ('O_iceareaN', iou, ntrec, a_noice, c1e4, c0)
      call putvars ('O_iceareaS', iou, ntrec, a_soice, c1e4, c0)
      call putvars ('O_snovolN', iou, ntrec, v_nsnow, c1e6, c0)
      call putvars ('O_snovolS', iou, ntrec, v_ssnow, c1e6, c0)
      call putvars ('A_satL', iou, ntrec, lsat, c1, c0)
      call putvars ('A_satO', iou, ntrec, osat, c1, c0)
      call putvars ('F_precipL', iou, ntrec, lprecip, p1, c0)
      call putvars ('F_precipO', iou, ntrec, oprecip, p1, c0)
      call putvars ('F_evapL', iou, ntrec, levap, p1, c0)
      call putvars ('F_evapO', iou, ntrec, oevap, p1, c0)
      call putvars ('F_solins', iou, ntrec, solins, c1e3, c0)
      call putvars ('F_upsens', iou, ntrec, upsens, c1e3, c0)
      call putvars ('F_uplwr', iou, ntrec, uplwr, c1e3, c0)
      call putvars ('F_outlwr', iou, ntrec, outlwr, c1e3, c0)
      call putvars ('F_dnswr', iou, ntrec, dnswr, c1e3, c0)
      call putvars ('F_absswr', iou, ntrec, absswr, c1e3, c0)
      call putvars ('F_netrad', iou, ntrec, netrad, c1e3, c0)
      call putvars ('A_albplt', iou, ntrec, palb, c1, c0)
      call putvars ('A_albatm', iou, ntrec, aalb, c1, c0)
      call putvars ('A_albsur', iou, ntrec, salb, c1, c0)
      call putvars ('A_albsurL', iou, ntrec, lsalb, c1, c0)
      call putvars ('A_albsurO', iou, ntrec, osalb, c1, c0)
      call putvars ('O_tempsur', iou, ntrec, sst, c1, c0)
      call putvars ('O_salsur', iou, ntrec, sss, p001, -p035)

      return
      end
