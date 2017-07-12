! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/fdift.h
!====================== include file "fdift.h" =========================

!     finite difference numerics for tracers
!=======================================================================

      T_i(i,k,j,n,ip) = t(i+ip,k,j,n,taum1)
      T_j(i,k,j,n,jp) = t(i,k,j+jp,n,taum1)
      dz_t2r(i,k,j) = dzt2r(k)
      dz_tr(i,k,j)  = dztr(k)
      dz_wtr(i,k,j) = dzwr(k)
      dx_t2r(i,k,j) = cstdxt2r(i,j)
      dx_tr(i,k,j)  = cstdxtr(i,j)
      dy_t2r(i,k,j) = cstdyt2r(jrow)
      dy_tr(i,k,j)  = cstdytr(jrow)

!-----------------------------------------------------------------------
!     advective terms
!-----------------------------------------------------------------------

      ADV_Tx(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*cstdxt2r(i,j)
      ADV_Ty(i,k,j,jrow,n) = (adv_f4n(i,k,j,n) - adv_f4n(i,k,j-1,n))
     &  *cstdyt2r(jrow)
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)

!     gent_mcwilliams isopycnal advective terms simulating the effect
!     of eddies on the isopycnals

      ADV_Txiso(i,k,j,n) = cstdxt2r(i,j)*(adv_vetiso(i,k,j)
     &  *(t(i+1,k,j,n,taum1) + t(i,k,j,n,taum1)) - adv_vetiso(i-1,k,j)
     &  *(t(i,k,j,n,taum1) + t(i-1,k,j,n,taum1)))
      ADV_Tyiso(i,k,j,jrow,n) = cstdyt2r(jrow)*(adv_vntiso(i,k,j)
     &  *(t(i,k,j+1,n,taum1) + t(i,k,j,n,taum1)) - adv_vntiso(i,k,j-1)
     &  *(t(i,k,j,n,taum1) + t(i,k,j-1,n,taum1)))
      ADV_Tziso(i,k,j) = dzt2r(k)*(adv_fbiso(i,k-1,j)-adv_fbiso(i,k,j))

!-----------------------------------------------------------------------
!     diffusive terms
!-----------------------------------------------------------------------

!     zonal component

      DIFF_Tx(i,k,j) = (diff_fe(i,  k,j)*tmask(i+1,k,j)
     &  - diff_fe(i-1,k,j)*tmask(i-1,k,j))*cstdxtr(i,j)

!     meridional component

      DIFF_Ty(i,k,j,jrow,n) = (diff_fn(i,k,j  )*tmask(i,k,j+1)
     &  - diff_fn(i,k,j-1)*tmask(i,k,j-1))*cstdytr(jrow)

!     vertical component

      DIFF_Tz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
     &  *(c1-aidif)
     &  + (diff_fbiso(i,k-1,j) - diff_fbiso(i,k,j))*dztr(k)
