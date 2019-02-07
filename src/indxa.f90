!*==INDXA.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK INDXA
SUBROUTINE INDXA(I,Ir,Idxa,Na)
  IMPLICIT NONE
  !*--INDXA5
  !*** Start of declarations inserted by SPAG
  REAL CNV, EPS
  INTEGER I, Idxa, IK, Ir, K, Na, NCMplx, NM, NPP
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  INDXA
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (INDXA-I)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  INDXA
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  INDXA
  Na = 2**Ir
  Idxa = I - Na + 1
  IF ( I>NM ) Na = 0
END SUBROUTINE INDXA
