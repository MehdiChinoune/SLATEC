!*==CSHCH.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CSHCH
SUBROUTINE CSHCH(Z,Csh,Cch)
  IMPLICIT NONE
  !*--CSHCH5
  !***BEGIN PROLOGUE  CSHCH
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESH and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CSHCH-A, ZSHCH-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
  !     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
  !
  !***SEE ALSO  CBESH, CBESK
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CSHCH
  COMPLEX Cch , Csh , Z
  REAL cchi , cchr , ch , cn , cshi , cshr , sh , sn , x , y
  !***FIRST EXECUTABLE STATEMENT  CSHCH
  x = REAL(Z)
  y = AIMAG(Z)
  sh = SINH(x)
  ch = COSH(x)
  sn = SIN(y)
  cn = COS(y)
  cshr = sh*cn
  cshi = ch*sn
  Csh = CMPLX(cshr,cshi)
  cchr = ch*cn
  cchi = sh*sn
  Cch = CMPLX(cchr,cchi)
END SUBROUTINE CSHCH
