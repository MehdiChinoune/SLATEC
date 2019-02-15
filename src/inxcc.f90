!DECK INXCC
SUBROUTINE INXCC(I,Ir,Idxc,Nc)
  IMPLICIT NONE
  REAL CNV, EPS
  INTEGER I, Idxc, IK, Ir, K, Nc, NCMplx, NM, NPP
  !***BEGIN PROLOGUE  INXCC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (INXCC-I)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CCBLK
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  INXCC
  COMMON /CCBLK / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  INXCC
  Nc = 2**Ir
  Idxc = I
  IF ( Idxc+Nc-1>NM ) Nc = 0
END SUBROUTINE INXCC
