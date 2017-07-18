! source file: /data/home/kai/dev/UVic2.9/nobio/updates/isopyc.h
!======================== include file "isopyc.h" ======================

!     isopycnal diffusion variables:

!     ahisop = isopycnal tracer mixing coefficient (cm**2/sec)
!     beta   = df/dy where f is related to coriolis force
!     drodx  = d(rho)/dx local to east face of T cell
!     drody  = d(rho)/dy local to north face of T cell
!     drodz  = d(rho)/dz local to bottom face of T cell
!     Ai_e   = diffusion coefficient on eastern face of T cell
!     Ai_n   = diffusion coefficient on northern face of T cell
!     Ai_bx  = diffusion coefficient on bottom face of T cell
!     Ai_by  = diffusion coefficient on bottom face of T cell

!     fisop  = structure function for isopycnal diffusion coefficient.
!     slmxr  = reciprocal of maximum allowable slope of isopycnals for
!              small angle approximation

      real alphai, betai, beta
      common /cisop_r/ alphai(imt,km,jmw), betai(imt,km,jmw), beta

      real addisop
      real ddxt, ddyt, ddzt, Ai_ez, Ai_nz, Ai_bx, Ai_by, K11, K22, K33
      real ahisop, fisop, slmxr, del_dm, s_dmr
      common /cisop_r/ ddxt(imt,km,jsmw:jemw,2)
      common /cisop_r/ ddyt(imt,km,1:jemw,2)
      common /cisop_r/ ddzt(imt,0:km,jmw,2)

      common /cisop_r/ Ai_ez(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ Ai_nz(imt,km,1:jemw,0:1,0:1)
      common /cisop_r/ Ai_bx(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ Ai_by(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ K11(imt,km,jsmw:jemw)
      common /cisop_r/ K22(imt,km,1:jemw)
      common /cisop_r/ K33(imt,km,jsmw:jemw)
      common /cisop_r/ ahisop, fisop(imt,jmt,km), slmxr
      common /cisop_r/ addisop(imt,km,jsmw:jemw)
      real delta_iso, s_minus, s_plus
      common /cisop_r/ delta_iso, s_minus, s_plus
!     adv_vetiso = zonal isopycnal mixing velocity computed at the
!                  center of the eastern face of the "t" cells
!     adv_vntiso = meridional isopycnal mixing velocity computed at
!                  the center of the northern face of the "t" cells
!                  (Note: this includes the cosine as in "adv_vnt")
!     adv_vbtiso = vertical isopycnal mixing velocity computed at the
!                  center of the top face of the "t" cells
!     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
!                  the bottom face of the "t" cells
!     athkdf = isopycnal thickness diffusivity (cm**2/sec)

      real athkdf, adv_vetiso, adv_vntiso, adv_vbtiso, adv_fbiso
      common /cisop_r/ athkdf
      common /cisop_r/ adv_vetiso(imt,km,jsmw:jemw)
      common /cisop_r/ adv_vntiso(imt,km,1:jemw)
      common /cisop_r/ adv_vbtiso(imt,0:km,jsmw:jemw)
      common /cisop_r/ adv_fbiso(imt,0:km,jsmw:jemw)

!     Define variables related to calculating K_gm mesoscale eddy
!     diffiusivity as outlined in Gent an McWilliams Paper (1989).
!     Further refinement from Eden 2009
!     *** NEED MORE COMMENTS HERE AND DETAIL ***
!     L_Rhi      = Rhines scale. Defined as sigma/beta. Where sigma is
!                  the Eady Growth rate of baroclinic instability
!     Lr         = 1st baroclinic Rossby Radius
!     c_eden     = Determined constant to ensure an average O_KGM =~ 800 m^2/s
!     kgm        = Isopycnal diffisivity constant
!     kgm_ave    = Average of Kgm. Mainly used to compute c_eden
!     kgm_sum
!     ahisop_var = Related to ahisop although this is a variable, vectorized
!                  version that changes with kgm and a constant
!     ahisop_sum
!     ahisop_ave
!     gridsum_area

      real drodxte, drodxbe
      real drodytn, drodybn
      real drodzte, drodzbe
      real drodztn, drodzbn

      integer countx, county
      real abs_grd_rho2, grd_rho_x, grd_rho_y, abs_drho_dz

      common /cisop_r/ drodxte(imt,km,jmt), drodxbe(imt,km,jmt)
      common /cisop_r/ drodytn(imt,km,jmt), drodybn(imt,km,jmt)
      common /cisop_r/ drodzte(imt,km,jmt), drodzbe(imt,km,jmt)
      common /cisop_r/ drodztn(imt,km,jmt), drodzbn(imt,km,jmt)

      !     Oleg and Geoff
      !     niso = number of indices in kgm. =2 for anisotropic GM coeff

      integer niso
      parameter (niso = 1)

      real Lm, Lr, L_Rhi, kgm, ahisop_var, gridsum_area
      real ahisop_sum, ahisop_ave, c_eden, coef, kgm_ave, kgm_sum, pii
      real stratif_int, clinic_int(niso), sum_zz
      real eddy_min, eddy_max

      common /kgm2d_r/ kgm(imt,jmt,niso), ahisop_var(imt,jmt,niso),
     &                 Lr(imt,jmt), L_Rhi(imt,jmt)

      real drodxe, drodze, drodyn, drodzn, drodxb, drodyb, drodzb
      real drodye, drodxn

!     statement functions

      drodxe(i,k,j,ip) =    alphai(i+ip,k,j)*ddxt(i,k,j,1) +
     &                      betai(i+ip,k,j)*ddxt(i,k,j,2)
      drodze(i,k,j,ip,kr) = alphai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,1) +
     &                      betai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,2)

      drodyn(i,k,j,jq) =    alphai(i,k,j+jq)*ddyt(i,k,j,1) +
     &                      betai(i,k,j+jq)*ddyt(i,k,j,2)
      drodzn(i,k,j,jq,kr) = alphai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,1) +
     &                      betai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,2)

      drodxb(i,k,j,ip,kr) = alphai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,1) +
     &                      betai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,2)
      drodyb(i,k,j,jq,kr) = alphai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,1) +
     &                      betai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,2)
      drodzb(i,k,j,kr) =    alphai(i,k+kr,j)*ddzt(i,k,j,1) +
     &                      betai(i,k+kr,j)*ddzt(i,k,j,2)

