!*==CDIV.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CDIV
SUBROUTINE CDIV(Ar,Ai,Br,Bi,Cr,Ci)
  IMPLICIT NONE
  !*--CDIV5
  !***BEGIN PROLOGUE  CDIV
  !***SUBSIDIARY
  !***PURPOSE  Compute the complex quotient of two complex numbers.
  !***LIBRARY   SLATEC
  !***TYPE      COMPLEX (CDIV-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
  !
  !***SEE ALSO  EISDOC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   811101  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CDIV
  REAL Ar , Ai , Br , Bi , Cr , Ci
  !
  REAL s , ars , ais , brs , bis
  !***FIRST EXECUTABLE STATEMENT  CDIV
  s = ABS(Br) + ABS(Bi)
  ars = Ar/s
  ais = Ai/s
  brs = Br/s
  bis = Bi/s
  s = brs**2 + bis**2
  Cr = (ars*brs+ais*bis)/s
  Ci = (ais*brs-ars*bis)/s
END SUBROUTINE CDIV
