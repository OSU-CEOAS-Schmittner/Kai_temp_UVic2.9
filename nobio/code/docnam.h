! source file: /usr/local/models/UVic_ESCM/2.9/source/mom/docnam.h
!====================== include file "docnam.h" ========================

!    info from docum.F that can be used elsewhere.
!    user specified tracer names are place into "trname" in docum

      character(12) :: trname
      common /docnam/ trname(nt)
