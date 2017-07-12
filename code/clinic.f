! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers8_npzd/updates/clinic.F
      subroutine clinic (joff, js, je, is, ie)

!=======================================================================
!     compute internal mode velocity components for rows js through je
!     in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer istrt, iend
      integer i, k, j, jrow, n, js, je, limit, joff, kb, is, ie

      real adv_ux, adv_uy, adv_uz, adv_metric, diff_ux, diff_uz, fx
      real diff_uy, diff_metric, coriolis, aprime, grav_rho0r, fxa
      real fxb, t1, t2, ambi_csur, del2, ambi_cst_dytr, detmr, unep

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "hmixc.h"
      include "emode.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"
      include "fdifm.h"

      parameter (istrt=2, iend=imt-1)

      real tempik(imt,km,jsmw:jmw)
      real baru(imt,jsmw:jemw,2)

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: these are currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,km
           do i=istrt-1,iend
              csudxur(i,j)  = csur(jrow)*dxur(i)
              csudxu2r(i,j) = csur(jrow)*dxur(i)*p5
              am_csudxtr(i,k,j) = visc_ceu(i,k,j)*csur(jrow)*dxtr(i+1)
           enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct the hydrostatic pressure gradients: 1 = dp/dx; 2 = dp/dy
!-----------------------------------------------------------------------

!     compute horizontal pressure gradient at the first level

      grav_rho0r = grav*rho0r
      do j=js,je
        jrow = j + joff
        fxa  = grav_rho0r*dzw(0)*csur(jrow)
        fxb  = grav_rho0r*dzw(0)*dyu2r(jrow)
        do i=istrt-1,iend
          t1              = rho(i+1,1,j+1) - rho(i  ,1,j)
          t2              = rho(i  ,1,j+1) - rho(i+1,1,j)
          grad_p(i,1,j,1) = (t1-t2)*fxa*dxu2r(i)
          grad_p(i,1,j,2) = (t1+t2)*fxb
        enddo
      enddo

!     compute the change in pressure gradient between levels

      do j=js,je+1
        do k=2,km
          do i=istrt-1,iend+1
            tempik(i,k,j) = rho(i,k-1,j) + rho(i,k,j)
          enddo
        enddo
      enddo

      do j=js,je
        jrow = j + joff
        fxa = grav_rho0r*csur(jrow)*p5
        fxb = grav_rho0r*dyu4r(jrow)
        do k=2,km
          do i=istrt-1,iend
            t1              = tempik(i+1,k,j+1) - tempik(i  ,k,j)
            t2              = tempik(i  ,k,j+1) - tempik(i+1,k,j)
            grad_p(i,k,j,1) = fxa*(t1-t2)*dzw(k-1)*dxu2r(i)
            grad_p(i,k,j,2) = fxb*(t1+t2)*dzw(k-1)
          enddo
        enddo
      enddo

!     integrate downward from the first level

      do j=js,je
        do k=1,kmm1
          do i=istrt-1,iend
            grad_p(i,k+1,j,1) = grad_p(i,k,j,1) + grad_p(i,k+1,j,1)
            grad_p(i,k+1,j,2) = grad_p(i,k,j,2) + grad_p(i,k+1,j,2)
          enddo
        enddo
      enddo

      do j=js,je
        call setbcx (grad_p(1,1,j,1), imt, km)
        call setbcx (grad_p(1,1,j,2), imt, km)
      enddo

!-----------------------------------------------------------------------
!     solve for one component of velocity at a time
!     n = 1 => zonal component
!     n = 2 => meridional component
!-----------------------------------------------------------------------

      do n=1,2

!-----------------------------------------------------------------------
!       calculate 2*advective flux (for speed) across east face of
!       "u" cells.
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              adv_fe(i,k,j) = adv_veu(i,k,j)*(u(i,  k,j,n,tau) +
     &                                        u(i+1,k,j,n,tau))
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       2*advective flux across northern face of "u" cells is built
!       into ADV_Uy. (It's done this way for performance issues)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       diffusive flux across east face of "u" cell
!       diffusive flux across north face of "u" cell
!-----------------------------------------------------------------------

!       build diffusive flux on eastern face of "u" cells

        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = am_csudxtr(i,k,j)*
     &                        (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1)
     &                        )
            enddo
          enddo
        enddo

!       diffusive flux on northern face of "u" cells is built
!       into DIFF_Uy

!-----------------------------------------------------------------------
!       calculate 2*advective flux (for speed) on bottom face of
!       "u" cell. also diffusive flux on bottom face of "u" cell
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,kmm1
            do i=istrt,iend
              adv_fb(i,k,j)  = adv_vbu(i,k,j)*(u(i,k,  j,n,tau) +
     &                                        u(i,k+1,j,n,tau))
              diff_fb(i,k,j) = visc_cbu(i,k,j)*dzwr(k)*
     &                         (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "u" cells for mixing
!       and advection.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            kb              = kmu(i,jrow)
            diff_fb(i,0,j)  = smf(i,j,n)
            diff_fb(i,kb,j) = bmf(i,j,n)
            adv_fb(i,0,j)   = adv_vbu(i,0,j)*(u(i,1,j,n,tau) +
     &                                        u(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbu(i,km,j)*u(i,km,j,n,tau)
          enddo
        enddo

!-----------------------------------------------------------------------
!       set source term for "u" cell
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            do i=istrt,iend
              source(i,k,j) = c0
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       solve for the internal mode part of du/dt at center of
!       "u" cells by neglecting the surface pressure gradients. use
!       statement functions to represent each component of the
!       calculation.
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) =
     &          (DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j,jrow,n) + DIFF_Uz(i,k,j)
     &          + DIFF_metric(i,k,j,jrow,n)
     &          - ADV_Ux(i,k,j) - ADV_Uy(i,k,j,jrow,n) - ADV_Uz(i,k,j)
     &          + ADV_metric(i,k,j,jrow,n)
     &          - grad_p(i,k,j,n) + CORIOLIS(i,k,j,jrow,n)
     &          + source(i,k,j)
     &          )*umask(i,k,j)
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with velocity component "n"
!-----------------------------------------------------------------------

        call diagc1 (joff, js, je, istrt, iend, n)

!-----------------------------------------------------------------------
!       construct the vertical average of du/dt for forcing
!       the barotropic equation
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = c0
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km
            fx = dzt(k)
            do i=istrt,iend
              zu(i,jrow,n) = zu(i,jrow,n) + u(i,k,j,n,taup1)*fx
            enddo
          enddo
        enddo

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = zu(i,jrow,n)*hr(i,jrow)
          enddo
        enddo

!-----------------------------------------------------------------------
!       end of velocity component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     compute "tau+1" velocities accounting for implicit part of the
!     coriolis term if treated implicitly. velocities are in error by an
!     arbitrary constant related to neglecting the unknown surface
!     pressure gradients
!-----------------------------------------------------------------------

      do n=1,2
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taum1)
     &                            + c2dtuv*u(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     subtract incorrect vertical means (related to ignoring horizontal
!     gradients of the surface pressure) to get pure internal modes.
!-----------------------------------------------------------------------

      do n=1,2
        do j=js,je
          do i=istrt,iend
            baru(i,j,n) = c0
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              baru(i,j,n) = baru(i,j,n) + u(i,k,j,n,taup1)*dzt(k)
            enddo
          enddo
        enddo
        do j=js,je
          jrow  = j + joff
          do i=istrt,iend
            baru(i,j,n) = baru(i,j,n)*hr(i,jrow)
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taup1)
     &                          - umask(i,k,j)*baru(i,j,n)
            enddo
          enddo
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics involving internal mode velocity at "tau+1"
!-----------------------------------------------------------------------

      call diagc2 (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!     filter velocity components at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filuv (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in clinic'
        stop '=>clinic'
      endif
      do j=js,je
        call setbcx (u(1,1,j,1,taup1), imt, km)
        call setbcx (u(1,1,j,2,taup1), imt, km)
      enddo

!-----------------------------------------------------------------------
!     if needed, construct the Ice S.B.C.(surface boundary conditions)
!     averaged over this segment
!-----------------------------------------------------------------------

      call isbcu (joff, js, je, istrt, iend, igu, igv)
      call asbcu (joff, js, je, istrt, iend, isu, isv)

      return
      end

      subroutine diagc1 (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     construct diagnostics which don`t require internal mode velocity
!     at "tau+1" for each velocity component "n"

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = (1,2) = (u,v) velocity component
!-----------------------------------------------------------------------

      implicit none

      integer n, j, js, je, jrow, joff, k, i, is, ie

      real dudx, ce, dudy, cn, dudz, cb, fx, weight

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "diag.h"
      include "diaga.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

      real temp(imt,km)

!-----------------------------------------------------------------------
!     diagnostic: accumulate global kinetic energy on "tau" velocity
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          fx = rho0*p5*csu(jrow)*dyu(jrow)
          do k=1,km
            do i=is,ie
              weight    = fx*dzt(k)*dxu(i)
              temp(i,k) = u(i,k,j,n,tau)**2*weight
            enddo
            do i=is,ie
              ektot(k,jrow) = ektot(k,jrow) + temp(i,k)
            enddo
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     diagnostic: integrate work done by the r.h.s. terms in the
!                  momentum equations.
!-----------------------------------------------------------------------

      if (glents .and. eots) call ge1 (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the momentum equations
!                 over specified regional volumes
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call utb1 (joff, js, je, is, ie, n)
      return
      end

      subroutine diagc2 (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!     construct diagnostics requiring internal mode velocity at "tau+1"
!     and those not dependent on velocity component fluxes.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!-----------------------------------------------------------------------

      implicit none

      integer joff, js, je, is, ie

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "scalar.h"
      include "switch.h"

!-----------------------------------------------------------------------
!     diagnostic: integrate work done by du/dt in the momentum equations
!                 the external mode part of "u" at "tau+1" will be
!                 accounted for after the external mode is solved.
!                 also, integrate the work done by buoyancy.
!-----------------------------------------------------------------------

      if (glents .and. eots) then
        call ge2 (joff, js, je, is, ie, kmt, kmu, c2dtuv, grav, rho0r)
      endif

!-----------------------------------------------------------------------
!     diagnostic: add du/dt and implicit coriolis terms to the integrals
!                 over specified volumes. the external mode parts will
!                 be accounted for after the external mode is solved.
!-----------------------------------------------------------------------

      if (trmbts .and. eots) then
        call utb2 (joff, js, je, is, ie, c2dtuv, acor)
      endif

      return
      end

      subroutine asbcu (joff, js, je, is, ie, iu, iv)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C.(surface boundary conditions)

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       iu   = index for u component
!       iv   = index for v component

!     reference: Pacanowski, R.C., Effect of Equatorial Currents
!                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
!-----------------------------------------------------------------------

      implicit none

      integer iu, iv, j, js, je, jrow, joff, i, is, ie

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

!     initialize S.B.C. at the beginning of each ocean segment
!     (do not alter values in land)

      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = c0
              sbc(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif

!     accumulate surface currents for the Atmos S.B.C. every time step

      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbc(i,jrow,iu) = sbc(i,jrow,iu) + p25*(
     &                       u(i,1,j,1,tau) + u(i-1,1,j,1,tau)
     &                     + u(i,1,j-1,1,tau) + u(i-1,1,j-1,1,tau))
            sbc(i,jrow,iv) = sbc(i,jrow,iv) + p25*(
     &                       u(i,1,j,2,tau) + u(i-1,1,j,2,tau)
     &                     + u(i,1,j-1,2,tau) + u(i-1,1,j-1,2,tau))
          enddo
        enddo
      endif

!     average the surface currents for the Atmos S.B.C. at the end
!     of each ocean segment. (do not alter values in land)

      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = rts*sbc(i,jrow,iu)
              sbc(i,jrow,iv) = rts*sbc(i,jrow,iv)
            endif
          enddo
        enddo
      endif

      return
      end

      subroutine isbcu (joff, js, je, is, ie, iu, iv)

!-----------------------------------------------------------------------
!     construct the Ice S.B.C.(surface boundary conditions)

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       iu   = index for u component
!       iv   = index for v component

!-----------------------------------------------------------------------

      implicit none

      integer iu, iv, j, js, je, jrow, joff, i, is, ie

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

!     initialize S.B.C. at the beginning of each ocean segment
!     (do not alter values in land)

      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = c0
              sbc(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif

!     accumulate geostrophic currents (level 2) for the Ice S.B.C.
!     every time step

      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbc(i,jrow,iu) = sbc(i,jrow,iu) + u(i,2,j,1,tau)
            sbc(i,jrow,iv) = sbc(i,jrow,iv) + u(i,2,j,2,tau)
          enddo
        enddo
      endif

!     average the currents for the Ice S.B.C. at the end of each ocean
!     segment. (do not alter values in land)

      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = rts*sbc(i,jrow,iu)
              sbc(i,jrow,iv) = rts*sbc(i,jrow,iv)
            endif
          enddo
        enddo
      endif

      return
      end

