!** DCBRT
REAL(8) FUNCTION DCBRT(X)
  !>
  !  Compute the cube root.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C2
  !***
  ! **Type:**      DOUBLE PRECISION (CBRT-S, DCBRT-D, CCBRT-C)
  !***
  ! **Keywords:**  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DCBRT(X) calculates the double precision cube root for
  ! double precision argument X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9PAK, D9UPAK

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  REAL(8) :: X
  INTEGER :: irem, iter, ixpnt, n
  REAL(8) :: y, cbrtsq
  REAL :: z
  REAL(8), PARAMETER :: cbrt2(5) = [ 0.62996052494743658238360530363911D0, &
    0.79370052598409973737585281963615D0, 1.0D0, &
    1.25992104989487316476721060727823D0, 1.58740105196819947475170563927231D0 ]
  INTEGER, PARAMETER :: niter = INT( 1.443*LOG(-.106*LOG(0.1*D1MACH(3))) )
  !* FIRST EXECUTABLE STATEMENT  DCBRT
  !
  DCBRT = 0.D0
  IF ( X==0.D0 ) RETURN
  !
  CALL D9UPAK(ABS(X),y,n)
  ixpnt = n/3
  irem = n - 3*ixpnt + 3
  !
  ! THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED
  ! TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF
  ! RELATIVE ERROR WITH 4.085 DIGITS ACCURACY.
  !
  z = REAL( y, 4 )
  DCBRT = .439581E0 + z*(.928549E0+z*(-.512653E0+z*.144586E0))
  !
  DO iter = 1, niter
    cbrtsq = DCBRT*DCBRT
    DCBRT = DCBRT + (y-DCBRT*cbrtsq)/(3.D0*cbrtsq)
  END DO
  !
  DCBRT = D9PAK(cbrt2(irem)*SIGN(DCBRT,X),ixpnt)
  !
END FUNCTION DCBRT
