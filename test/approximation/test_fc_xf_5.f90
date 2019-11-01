MODULE TEST_FC_XF_5
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE FCQX()
    !> Quick check for FC.
    !***
    ! **Description:**
    !
    !   Quick check subprogram for the subroutine FC.
    !
    !   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
    !   its first two derivatives, and probable error curve.
    !
    !   Use subprogram FC to obtain the constrained cubic B-spline
    !   representation of the curve.
    !
    !   The values of the coefficients of the B-spline as computed by FC
    !   and the values of the fitted curve as computed by BVALU in the
    !   de Boor package are tested for accuracy with the expected values.
    !   See the example program in the report sand78-1291, pp. 22-27.
    !
    !   The dimensions in the following arrays are as small as possible for
    !   the problem being solved.

    USE approximation, ONLY : FC
    !     .. Local Scalars ..
    REAL(SP) :: t
    INTEGER :: i, idigit, l, mode, nconst, ndeg
    !     .. Local Arrays ..
    REAL(SP) :: coeff(9), w(529), xconst(11), yconst(11)
    INTEGER :: iw(30), nderiv(11)
    !     .. Data statements ..
    !
    INTEGER, PARAMETER :: ndata = 9, nord = 4, nbkpt = 13, last = 10
    REAL(SP), PARAMETER :: xdata(9) = [ 0.15_SP, 0.27_SP, 0.33_SP, 0.40_SP, 0.43_SP, &
      0.47_SP, 0.53_SP, 0.58_SP, 0.63_SP ]
    REAL(SP), PARAMETER :: ydata(9) = [ 0.025_SP, 0.05_SP, 0.13_SP, 0.27_SP, 0.37_SP, &
      0.47_SP, 0.64_SP, 0.77_SP, 0.87_SP ]
    REAL(SP), PARAMETER :: sddata(9) = 0.015_SP
    REAL(SP), PARAMETER :: bkpt(13) = [ -0.6_SP, -0.4_SP, -0.2_SP, 0._SP, 0.2_SP, &
      0.4_SP, 0.6_SP, 0.8_SP, 0.9_SP, 1._SP, 1.1_SP, 1.2_SP, 1.3_SP ]
    !* FIRST EXECUTABLE STATEMENT  FCQX
    ndeg = nord - 1
    !
    !     Write the various constraints for the fitted curve.
    !
    nconst = 0
    t = bkpt(nord)
    !
    !     Constrain function to be zero at left-most breakpoint.
    !
    nconst = nconst + 1
    xconst(nconst) = t
    yconst(nconst) = 0._SP
    nderiv(nconst) = 2 + 4*0
    !
    !     Constrain first derivative to be nonnegative at left-most
    !     breakpoint.
    !
    nconst = nconst + 1
    xconst(nconst) = t
    yconst(nconst) = 0._SP
    nderiv(nconst) = 1 + 4*1
    !
    !  Constrain second derivatives to be nonnegative at left set of breakpoints.
    !
    DO i = 1, 3
      l = ndeg + i
      t = bkpt(l)
      nconst = nconst + 1
      xconst(nconst) = t
      yconst(nconst) = 0._SP
      nderiv(nconst) = 1 + 4*2
    END DO
    !
    !     Constrain function value at right-most breakpoint to be one.
    !
    nconst = nconst + 1
    t = bkpt(last)
    xconst(nconst) = t
    yconst(nconst) = 1._SP
    nderiv(nconst) = 2 + 4*0
    !
    !     Constrain slope to agree at left- and right-most breakpoints.
    !
    nconst = nconst + 1
    xconst(nconst) = bkpt(nord)
    yconst(nconst) = bkpt(last)
    nderiv(nconst) = 3 + 4*1
    !
    !     Constrain second derivatives to be nonpositive at right set of
    !     breakpoints.
    !
    DO i = 1, 4
      nconst = nconst + 1
      l = last - 4 + i
      xconst(nconst) = bkpt(l)
      yconst(nconst) = 0._SP
      nderiv(nconst) = 0 + 4*2
    END DO
    !
    idigit = -4
    !
    !     Declare amount of working storage allocated to FC.
    !
    iw(1) = 529
    iw(2) = 30
    !
    !  Set mode to indicate a new problem and request the variance function.
    !
    mode = 2
    !
    !     Trigger error conditions.
    !
    iw(1) = 10
    CALL FC(ndata,xdata,ydata,sddata,nord,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    !
    RETURN
  END SUBROUTINE FCQX
  !
END MODULE TEST_FC_XF_5
!
PROGRAM MAIN
  USE TEST_FC_XF_5
  IMPLICIT NONE
  !
  CALL FCQX()
  !
END PROGRAM MAIN