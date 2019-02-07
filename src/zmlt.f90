!*==ZMLT.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZMLT
SUBROUTINE ZMLT(Ar,Ai,Br,Bi,Cr,Ci)
  IMPLICIT NONE
  !*--ZMLT5
  !***BEGIN PROLOGUE  ZMLT
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (ZMLT-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     DOUBLE PRECISION COMPLEX MULTIPLY, C=A*B.
  !
  !***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZMLT
  REAL(8) :: Ar , Ai , Br , Bi , Cr , Ci , ca , cb
  !***FIRST EXECUTABLE STATEMENT  ZMLT
  ca = Ar*Br - Ai*Bi
  cb = Ar*Bi + Ai*Br
  Cr = ca
  Ci = cb
END SUBROUTINE ZMLT
