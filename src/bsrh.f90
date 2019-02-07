!*==BSRH.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK BSRH
FUNCTION BSRH(Xll,Xrr,Iz,C,A,Bh,F,Sgn)
  IMPLICIT NONE
  !*--BSRH5
  !*** Start of declarations inserted by SPAG
  REAL A, Bh, BSRH, C, CNV, dx, EPS, F, Sgn, x, xl, Xll, xr, &
    Xrr
  INTEGER IK, Iz, K, NCMplx, NM, NPP
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BSRH
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  BSRH
  DIMENSION A(*), C(*), Bh(*)
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  BSRH
  xl = Xll
  xr = Xrr
  dx = .5*ABS(xr-xl)
  100  x = .5*(xl+xr)
  IF ( Sgn*F(x,Iz,C,A,Bh)<0 ) THEN
    xl = x
  ELSEIF ( Sgn*F(x,Iz,C,A,Bh)==0 ) THEN
    BSRH = .5*(xl+xr)
    GOTO 99999
  ELSE
    xr = x
  ENDIF
  dx = .5*dx
  IF ( dx>CNV ) GOTO 100
  BSRH = .5*(xl+xr)
  99999 CONTINUE
  END FUNCTION BSRH
