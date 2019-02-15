!DECK ZEXP
SUBROUTINE ZEXP(Ar,Ai,Br,Bi)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ZEXP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and
  !            ZBIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (ZEXP-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     DOUBLE PRECISION COMPLEX EXPONENTIAL FUNCTION B=EXP(A)
  !
  !***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZEXP
  REAL(8) :: Ar, Ai, Br, Bi, zm, ca, cb
  !***FIRST EXECUTABLE STATEMENT  ZEXP
  zm = EXP(Ar)
  ca = zm*COS(Ai)
  cb = zm*SIN(Ai)
  Br = ca
  Bi = cb
END SUBROUTINE ZEXP
