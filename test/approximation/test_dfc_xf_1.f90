MODULE TEST_DFC_XF_1
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DFCQX()
    !> Quick check for DFC.
    !***
    ! **Description:**
    !
    !   Quick check subprogram for the subroutine DFC.
    !
    !   Fit discrete data by an S-shaped curve.  Evaluate the fitted curve,
    !   its first two derivatives, and probable error curve.
    !
    !   Use subprogram DFC to obtain the constrained cubic B-spline
    !   representation of the curve.
    !
    !   The values of the coefficients of the B-spline as computed by DFC
    !   and the values of the fitted curve as computed by DBVALU in the
    !   de Boor package are tested for accuracy with the expected values.
    !   See the example program in the report sand78-1291, pp. 22-27.
    !
    !   The dimensions in the following arrays are as small as possible for
    !   the problem being solved.

    USE approximation, ONLY : DFC
    !     .. Local Scalars ..
    REAL(DP) :: t
    INTEGER :: i, idigit, l, mode, nconst, ndeg
    !     .. Local Arrays ..
    REAL(DP) :: coeff(9), w(529), xconst(11), yconst(11)
    INTEGER :: iw(30), nderiv(11)
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, REAL, SQRT
    !     .. Data statements ..
    !
    INTEGER, PARAMETER :: ndata = 9, nord = 4, nbkpt = 13, last = 10
    REAL(DP), PARAMETER :: xdata(9) = [ 0.15_DP, 0.27_DP, 0.33_DP, 0.40_DP, 0.43_DP, &
      0.47_DP, 0.53_DP, 0.58_DP, 0.63_DP ]
    REAL(DP), PARAMETER :: ydata(9) = [ 0.025_DP, 0.05_DP, 0.13_DP, 0.27_DP, 0.37_DP, &
      0.47_DP, 0.64_DP, 0.77_DP, 0.87_DP ]
    REAL(DP), PARAMETER :: sddata(9)  = 0.015_DP
    REAL(DP), PARAMETER :: bkpt(13) = [ -0.6_DP, -0.4_DP, -0.2_DP, 0._DP, 0.2_DP, &
      0.4_DP, 0.6_DP, 0.8_DP, 0.9_DP, 1._DP, 1.1_DP, 1.2_DP, 1.3_DP ]
    !* FIRST EXECUTABLE STATEMENT  DFCQX
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
    yconst(nconst) = 0._DP
    nderiv(nconst) = 2 + 4*0
    !
    !  Constrain first derivative to be nonnegative at left-most breakpoint.
    !
    nconst = nconst + 1
    xconst(nconst) = t
    yconst(nconst) = 0._DP
    nderiv(nconst) = 1 + 4*1
    !
    !  Constrain second derivatives to be nonnegative at left set of breakpoints.
    !
    DO i = 1, 3
      l = ndeg + i
      t = bkpt(l)
      nconst = nconst + 1
      xconst(nconst) = t
      yconst(nconst) = 0._DP
      nderiv(nconst) = 1 + 4*2
    END DO
    !
    !     Constrain function value at right-most breakpoint to be one.
    !
    nconst = nconst + 1
    t = bkpt(last)
    xconst(nconst) = t
    yconst(nconst) = 1._DP
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
      yconst(nconst) = 0._DP
      nderiv(nconst) = 0 + 4*2
    END DO
    !
    idigit = -4
    !
    !     Declare amount of working storage allocated to DFC.
    !
    iw(1) = 529
    iw(2) = 30
    !
    !  Set mode to indicate a new problem and request the variance function.
    !
    mode = 2
    !
    !  Trigger error conditions.
    !
    CALL DFC(ndata,xdata,ydata,sddata,0,nbkpt,bkpt,nconst,xconst,yconst,&
      nderiv,mode,coeff,w,iw)
    !
    RETURN
  END SUBROUTINE DFCQX
  !
END MODULE TEST_DFC_XF_1
!
PROGRAM MAIN
  USE TEST_DFC_XF_1
  IMPLICIT NONE
  !
  CALL DFCQX()
  !
END PROGRAM MAIN