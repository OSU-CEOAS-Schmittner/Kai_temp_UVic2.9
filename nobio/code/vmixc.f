! source file: /data/home/kai/dev/UVic2.9/nobio/updates/vmixc.F
      subroutine vmixc (joff, js, je, is, ie)

!=======================================================================
!     set viscosity coefficient on bottom face of "u" cells
!     set diffusion coefficient on bottom face of "t" cells

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, istrt, is, iend, ie, jstrt
      integer jend, jrow, joff, ks, k1

      real zn2, hab, zkappa, qk1, qo1, edr, q2

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"
      include "isopyc.h"
      include "tidal_kv.h"
      include "diag.h"
      include "grdvar.h"
      include "levind.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude and latitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)
      jstrt = max(2,js-1)
      jend  = je-1

!-----------------------------------------------------------------------
!     constant vertical mixing coefficients
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = j + joff
        do i=istrt,iend
          if (abs(tlat(i,jrow)) .lt. 30.) then
	    qk1 = 0.33
	    qo1 = 0.33
	  else
	    qk1 = 1.
	    qo1 = 1.
	  endif
          if (abs(tlat(i,jrow)) .lt. 70.) then
            q2 = 0.33
	  else
            q2 = 1.
	  endif
c          q2=0.5
          do k=1,kmt(i,jrow)-1
            visc_cbu(i,k,j) = kappa_m

!           calculate N^2 = -g/rho drhodz on bottom of cell face
!           (where K33 and diff_cbt = kappa_h are defined). Note that
!           N2 is not guaranteed to be positive. If instability occurs,
!           convective adjustment will eliminate it.
!           drodzb is defined in isopyc.h

!           ZN2 is defined on T-cell bottom (zw pt)
            ZN2 = max(-gravrho0r*drodzb(i,k,j,0),1e-8)

!           subgrid-scale scheme: sum over all levels below k
            edr = 0.
            do k1=k+1,kmt(i,jrow)
!              height above bottom
               hab = zw(k) - zw(k1)
               edr = edr + (q2*(edrm2(i,k1,jrow)+edrs2(i,k1,jrow))
     &                   + qk1*edrk1(i,k1,jrow) + qo1*edro1(i,k1,jrow))
     &                   *exp(hab*zetar)/(1-exp(-zetar*zw(k1)))
            enddo
            zkappa = ogamma*edr/ZN2

! Andreas: high mixing in Southern Ocean below 500 m as in observations
c            zkappa = max(zkappa,tanh((zw(k)-50000.)/10000)*
c     &                          (1-tanh((tlat(i,j+joff)+40)/8))/2)
!           limit diff_cbt
            diff_cbt(i,k,j) = max(kappa_h, min(100., zkappa + kappa_h))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     Add K33 component to vertical diffusion coefficient
!-----------------------------------------------------------------------

      do j=jstrt,jend
        do i=istrt,iend
          do k=1,km
            diff_cbt(i,k,j) = diff_cbt(i,k,j) + K33(i,k,j)
          enddo
        enddo
      enddo

      return
      end
