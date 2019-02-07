!*==ZSHCH.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZSHCH
SUBROUTINE ZSHCH(Zr,Zi,Cshr,Cshi,Cchr,Cchi)
  IMPLICIT NONE
  !*--ZSHCH5
  !***BEGIN PROLOGUE  ZSHCH
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CSHCH-A, ZSHCH-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
  !     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
  !
  !***SEE ALSO  ZBESH, ZBESK
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZSHCH
  !
  DOUBLE PRECISION Cchi , Cchr , ch , cn , Cshi , Cshr , sh , sn , Zi , Zr
  !***FIRST EXECUTABLE STATEMENT  ZSHCH
  sh = SINH(Zr)
  ch = COSH(Zr)
  sn = SIN(Zi)
  cn = COS(Zi)
  Cshr = sh*cn
  Cshi = ch*sn
  Cchr = ch*cn
  Cchi = sh*sn
END SUBROUTINE ZSHCH
