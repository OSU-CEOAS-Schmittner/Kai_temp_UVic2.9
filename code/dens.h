! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/dens.h
!====================== include file "dens.h" ==========================

!-----------------------------------------------------------------------
!     statement function
!-----------------------------------------------------------------------

      dens (tq, sq, k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &                   (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &                   (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq

      drodt (tq, sq, k) = c(k,1) + (c(k,4) + c(k,7)*sq)*sq + (2.0*c(k,3)
     &                  + 2.0*c(k,8)*sq + 3.0*c(k,6)*tq)*tq

      drods (tq, sq, k) = (c(k,4) + 2.0*c(k,7)*sq + c(k,8)*tq)*tq
     &                  + c(k,2) + (2.0*c(k,5) + 3.0*c(k,9)*sq)*sq
