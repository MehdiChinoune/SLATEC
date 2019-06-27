!** CBRT
REAL(SP) ELEMENTAL FUNCTION CBRT(X)
  !> Compute the cube root.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C2
  !***
  ! **Type:**      SINGLE PRECISION (CBRT-S, DCBRT-D, CCBRT-C)
  !***
  ! **Keywords:**  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CBRT(X) calculates the cube root of X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, R9PAK, R9UPAK

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  INTEGER :: irem, iter, ixpnt, n
  REAL(SP) :: cbrtsq, y
  REAL(SP), PARAMETER :: cbrt2(5) = [ 0.62996052494743658_SP,  0.79370052598409974_SP, &
    1._SP, 1.25992104989487316_SP, 1.58740105196819947_SP ]
  INTEGER, PARAMETER :: niter = INT( 1.443*LOG(-.106*LOG(0.1_SP*R1MACH(3))) ) + 1
  !* FIRST EXECUTABLE STATEMENT  CBRT
  !
  IF( X==0._SP ) THEN
    CBRT = 0._SP
    RETURN
  END IF
  !
  CALL R9UPAK(ABS(X),y,n)
  ixpnt = n/3
  irem = n - 3*ixpnt + 3
  !
  ! THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED
  ! TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF
  ! RELATIVE ERROR WITH 4.085 DIGITS ACCURACY.
  !
  CBRT = .439581_SP + y*(.928549_SP+y*(-.512653_SP+y*.144586E0_SP))
  !
  DO iter = 1, niter
    cbrtsq = CBRT*CBRT
    CBRT = CBRT + (y-CBRT*cbrtsq)/(3._SP*cbrtsq)
  END DO
  !
  CBRT = R9PAK(cbrt2(irem)*SIGN(CBRT,X),ixpnt)
  !
END FUNCTION CBRT