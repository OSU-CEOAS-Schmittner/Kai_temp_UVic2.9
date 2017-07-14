! source file: /usr/local/models/UVic_ESCM/2.9/source/common/emode.h
!======================= include file "emode.h" ========================

!     variables for  external mode

!     psi   = stream function (,,1) is for tau; (,,2) is for tau-1
!     zu    = vertically averaged forcing from momentum equations
!             (,,1) is zonal and (,,2) is meridional component
!     ztd   = curl of "zu" for the stream function equation
!     ptd   = time change of stream function
!     h     = depth over "u" points
!     hr    = reciprocal depth over "u" points
!     res   = residual from elliptic solver

!     map   = land mass map distinguishing, ocean, land, and perimeters

!     mxscan  = max number of allowable scans for Poisson solvers
!     mscan   = actual number of scans taken by Poisson solvers
!     sor     = successive over-relaxation constant
!     tolrsf  = tolerance for stream function calculation.
!               the solution is halted when it is within "tolrsf"
!               of the "true" solution assuming geometric convergence.
!     tolrsp  = tolerance for surface pressure calculation
!               the solution is halted when it is within "tolrsp"
!               of the "true" solution assuming geometric convergence.
!     tolrfs  = tolerance for implicit free surface calculation
!               the solution is halted when it is within "tolrfs"
!               of the "true" solution assuming geometric convergence.
!     esterr  = estimated maximum error in elliptic solver assuming
!               geometric convergence

!     nisle = number of land masses
!     nippts= number of land mass perimeter points
!     iperm = "i" coordinate for the land mass perimeter point
!     jperm = "j" coordinate for the land mass perimeter point
!     iofs  = offset for indexing into the land mass perimeter points
!     imask = controls whether calculations get done on perimeters
!     set mask for land mass perimeters on which to perform calculations
!     imask(-n) = .false.  [no equations ever on dry land mass n]
!     imask(0)  = .true.   [equations at all mid ocean points]
!     imask(n)  = .true./.false [controls whether there will be
!                                equations on the ocean perimeter of
!                                land mass n]
!     note: land mass 1 is the northwest-most land mass
!     for the numbering of the other landmasses, see generated map(i,j)

      character(16) :: variable
      common /emode_c/ variable

      logical converged, imask
      common /emode_l/ converged, imask (-mnisle:mnisle)

      integer mxscan, mscan, map, nippts, iofs, iperm
      integer jperm, nisle, imain

      common /emode_i/ mxscan, mscan
      common /emode_i/ map(imt,jmt)
      common /emode_i/ nippts(mnisle), iofs(mnisle), iperm(maxipp)
      common /emode_i/ jperm(maxipp), nisle, imain

      real tolrsf, tolrsp, tolrfs, sor, esterr, ptd, res, hr
      real h, zu, psi, ztd, alph, gam, theta, apgr, uhat, pguess
      real ps, divf, ubar, ubarm1

      common /emode_r/ tolrsf, tolrsp, tolrfs, sor, esterr
      common /emode_r/ ptd(imt,jmt), res(imt,jmt), hr(imt,jmt)
      common /emode_r/ h(imt,jmt), zu(imt,jmt,2)
      common /emode_r/ psi(imt,jmt,2), ztd(imt,jmt)
