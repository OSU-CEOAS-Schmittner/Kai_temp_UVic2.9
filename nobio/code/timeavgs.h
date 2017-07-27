! source file: /data/home/kai/dev/UVic2.9/nobio/updates/timeavgs.h
! source file: /Users/dkeller/Desktop/UVic_ESCM/2.9/source/mom/timeavgs.h
!====================== include file "timeavgs.h" ======================

!     imtav =  # of longitudes used for the time averages grid
!     jmtav =  # of latitudes used for the time averages grid
!     kmav  =  # of levels used for the time averages grid

      integer imtav, jmtav, kmav
      parameter (imtav=imt, jmtav=jmt-2, kmav=km)
      real ta_vetiso, ta_vntiso, ta_vbtiso
      common /ta_gm_r/ ta_vetiso(imt,km,jmt), ta_vntiso(imt,km,jmt)
      common /ta_gm_r/ ta_vbtiso(imt,km,jmt)

      real ta_kgm
      common /ta_gm_r/ ta_kgm(imt,km,jmt,1)

      integer nta_conv
      common /ta_conv_i/ nta_conv

      real ta_totalk, ta_vdepth, ta_pe
      common /ta_conv_r/ ta_totalk(imt,jmt), ta_vdepth(imt,jmt)
      common /ta_conv_r/ ta_pe(imt,jmt)
      integer nta_sscar
      common /ta_car_i/ nta_sscar

      real ta_rfeorgads, ta_rfecol
      real ta_rchl, ta_rchl_D
      real ta_rdeffe, ta_rremife, ta_rexpofe, ta_rfeprime
      real ta_rfesed, ta_rbfe
      real ta_rdeffe_C

      common /ta_npzd_r/ ta_rfeorgads(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfecol(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeffe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremife(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpofe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfeprime(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfesed(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rbfe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeffe_C(imt,kpzd,jmt)
