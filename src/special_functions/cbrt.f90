!** CBRT
REAL FUNCTION CBRT(X)
  !>
  !***
  !  Compute the cube root.
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

  REAL cbrtsq, X, y
  INTEGER irem, iter, ixpnt, n
  REAL, PARAMETER :: cbrt2(5) = [ 0.62996052494743658E0,  0.79370052598409974E0, &
    1.0E0, 1.25992104989487316E0, 1.58740105196819947E0 ]
  INTEGER :: niter = 0
  !* FIRST EXECUTABLE STATEMENT  CBRT
  IF ( niter==0 ) niter = INT( 1.443*LOG(-.106*LOG(0.1*R1MACH(3))) ) + 1
  !
  CBRT = 0.0
  IF ( X==0. ) RETURN
  !
  CALL R9UPAK(ABS(X),y,n)
  ixpnt = n/3
  irem = n - 3*ixpnt + 3
  !
  ! THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED
  ! TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF
  ! RELATIVE ERROR WITH 4.085 DIGITS ACCURACY.
  !
  CBRT = .439581E0 + y*(.928549E0+y*(-.512653E0+y*.144586E0))
  !
  DO iter = 1, niter
    cbrtsq = CBRT*CBRT
    CBRT = CBRT + (y-CBRT*cbrtsq)/(3.0*cbrtsq)
  END DO
  !
  CBRT = R9PAK(cbrt2(irem)*SIGN(CBRT,X),ixpnt)
  !
END FUNCTION CBRT
