!DECK BCRH
REAL FUNCTION BCRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BCRH
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CCBLK
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  BCRH
  REAL A, Bh, C, CNV, dx, EPS, F, Sgn, x, xl, Xll, xr, Xrr
  INTEGER IK, Iz, K, NCMplx, NM, NPP
  DIMENSION A(*), C(*), Bh(*)
  COMMON /CCBLK / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  BCRH
  xl = Xll
  xr = Xrr
  dx = .5*ABS(xr-xl)
  100  x = .5*(xl+xr)
  IF ( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF ( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BCRH = .5*(xl+xr)
    RETURN
  ELSE
    xr = x
  ENDIF
  dx = .5*dx
  IF ( dx>CNV ) GOTO 100
  BCRH = .5*(xl+xr)
  RETURN
END FUNCTION BCRH
