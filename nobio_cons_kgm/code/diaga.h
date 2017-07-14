! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/diaga.h
!====================== include file "diaga.h" =========================

!     variables used for computing diagnostics:

!     totalk   = total number of levels involved in convection
!     vdepth   = ventilation depth (cm)
!     pe       = potential energy lost due to explicit convection (g/s2)

      real totalk, vdepth, pe
      common /cdiaga_r/ totalk(imt,jsmw:jemw), vdepth(imt,jsmw:jemw)
      common /cdiaga_r/ pe(imt,jsmw:jemw)
