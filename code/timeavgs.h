! source file: /raid24/aschmitt/UVic2.9/karin/mwc15_npzd_fe_n_c13_alk_caco3/updates/timeavgs.h
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
      real ta_rnpp, ta_rgraz, ta_rmorp, ta_rmorpt, ta_rmorz, ta_rexcr
      real ta_rremi, ta_rexpo, ta_rnpp_D, ta_rgraz_D, ta_rmorpt_D
      real ta_rnfix, ta_rwcdeni, ta_rexpocal, ta_rprocal, ta_rgraz_Z
      real ta_rgraz_Det, ta_ravej, ta_ravej_D, ta_rgmax, ta_rno3P
      real ta_rpo4P, ta_rpo4_D, ta_kpipe, ta_rsedrr, ta_rbdeni
      real ta_rmorp_D, ta_rnpp_dop, ta_rnpp_D_dop, ta_rnpp_C_dop
      real ta_rexpocaco3, ta_rdissl, ta_rimpocaco3, ta_rcalatt
      real ta_rmorpt_C, ta_rnpp_C, ta_rgraz_C, ta_rmorp_C

      common /ta_npzd_r/ ta_rnpp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpocaco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rimpocaco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rdissl(imt,km,jmt)
      common /ta_npzd_r/ ta_rcalatt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rprocal(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp_dop(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp_C_dop(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexcr(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremi(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpocal(imt,km,jmt)
      common /ta_npzd_r/ ta_rsedrr(imt,jmt)
      common /ta_npzd_r/ ta_rnpp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp_D_dop(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnfix(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rwcdeni(imt,km,jmt)
      common /ta_npzd_r/ ta_rbdeni(imt,km,jmt)
      common /ta_npzd_r/ ta_rgraz_Det(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_Z(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ravej(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ravej_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgmax(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rno3P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_kpipe(imt,jmt)
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
