!** DCBRT
REAL(8) FUNCTION DCBRT(X)
  IMPLICIT NONE
  !>
  !***
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
  
  INTEGER irem, iter, ixpnt, n, niter
  REAL z
  REAL(8) :: X, cbrt2(5), y, cbrtsq, D9PAK, D1MACH
  SAVE cbrt2, niter
  DATA cbrt2(1)/0.62996052494743658238360530363911D0/
  DATA cbrt2(2)/0.79370052598409973737585281963615D0/
  DATA cbrt2(3)/1.0D0/
  DATA cbrt2(4)/1.25992104989487316476721060727823D0/
  DATA cbrt2(5)/1.58740105196819947475170563927231D0/
  DATA niter/0/
  !* FIRST EXECUTABLE STATEMENT  DCBRT
  IF ( niter==0 ) niter = INT( 1.443*LOG(-.106*LOG(0.1*REAL(D1MACH(3)))) )
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
  ENDDO
  !
  DCBRT = D9PAK(cbrt2(irem)*SIGN(DCBRT,X),ixpnt)
  !
END FUNCTION DCBRT
