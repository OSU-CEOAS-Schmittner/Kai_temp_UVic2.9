! source file: /raid24/aschmitt/UVic2.9/karin/mobi_with_calcifiers7_nobio/updates/fdifm.h
!====================== include file "fdifm.h" =========================

!     finite difference numerics for momentum
!=======================================================================

!-----------------------------------------------------------------------
!     advective terms
!-----------------------------------------------------------------------

      ADV_Ux(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*csudxu2r(i,j)
      ADV_Uy(i,k,j,jrow,n) = (adv_vnu(i,k,j)*(u(i,k,j,n,tau)
     &  + u(i,k,j+1,n,tau)) - adv_vnu(i,k,j-1)*(u(i,k,j-1,n,tau)
     &  + u(i,k,j,n,tau)))*csudyu2r(jrow)
      ADV_Uz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
      ADV_metric(i,k,j,jrow,n) = advmet(jrow,n)*u(i,k,j,1,tau)
     &  *u(i,k,j,3-n,tau)

!-----------------------------------------------------------------------
!     diffusive terms
!-----------------------------------------------------------------------

      DIFF_Ux(i,k,j) = (diff_fe(i,k,j) - diff_fe(i-1,k,j))
     &  *csudxur(i,j)
      DIFF_Uz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
      DIFF_Uy(i,k,j,jrow,n) = amc_north(i,k,jrow)*(u(i,k,j+1,n,taum1)
     &  - u(i,k,j,n,taum1)) - amc_south(i,k,jrow)*(u(i,k,j,n,taum1)
     &  - u(i,k,j-1,n,taum1))

!-----------------------------------------------------------------------
!     metric term
!-----------------------------------------------------------------------

      DIFF_metric(i,k,j,jrow,n) = am3(jrow)*u(i,k,j,n,taum1)
     &  + am4(jrow,n)*dxmetr(i)*(u(i+1,k,j,3-n,taum1)
     &  - u(i-1,k,j,3-n,taum1))

!-----------------------------------------------------------------------
!     coriolis term
!-----------------------------------------------------------------------

      CORIOLIS(i,k,j,jrow,n) = cori(i,jrow,n)*u(i,k,j,3-n,tau)
