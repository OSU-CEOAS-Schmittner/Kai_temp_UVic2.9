! source file: /usr/local/models/UVic_ESCM/2.9/source/embm/embmbc.F
      subroutine embmbc (data)

!=======================================================================
!     set boundary values

!     input:
!       data = array to be set
!=======================================================================

      implicit none

      integer i, j

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"

      real data(imt,jmt)

!     set cyclic east west

      do j=1,jmt
        data(1,j) = data(imtm1,j)
        data(imt,j) = data(2,j)
      enddo

!     set closed North South

      do i=1,imt
        data(i,1) = data(i,2)
        data(i,jmt) = data(i,jmtm1)
      enddo

      return
      end
