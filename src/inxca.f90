!DECK INXCA
SUBROUTINE INXCA(I,Ir,Idxa,Na)
  IMPLICIT NONE
  REAL CNV, EPS
  INTEGER I, Idxa, IK, Ir, K, Na, NCMplx, NM, NPP
  !***BEGIN PROLOGUE  INXCA
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (INXCA-I)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CCBLK
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  INXCA
  COMMON /CCBLK / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  INXCA
  Na = 2**Ir
  Idxa = I - Na + 1
  IF ( I>NM ) Na = 0
END SUBROUTINE INXCA
