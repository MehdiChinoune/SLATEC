!DECK INDXB
SUBROUTINE INDXB(I,Ir,Idx,Idp)
  IMPLICIT NONE
  REAL CNV, EPS
  INTEGER I, id, Idp, Idx, IK, ipl, Ir, izh, K, NCMplx, NM, NPP
  !***BEGIN PROLOGUE  INDXB
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      INTEGER (INDXB-I)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    CBLKT
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920422  Added statement so IDX would always be defined.  (WRB)
  !***END PROLOGUE  INDXB
  !
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !***FIRST EXECUTABLE STATEMENT  INDXB
  Idx = max(I,1)
  Idp = 0
  IF ( Ir<0 ) RETURN
  IF ( Ir==0 ) THEN
    IF ( I>NM ) RETURN
    Idx = I
    Idp = 1
    RETURN
  ELSE
    izh = 2**Ir
    id = I - izh - izh
    Idx = id + id + (Ir-1)*IK + Ir + (IK-I)/izh + 4
    ipl = izh - 1
    Idp = izh + izh - 1
    IF ( I-ipl<=NM ) THEN
      IF ( I+ipl>NM ) Idp = NM + ipl - I + 1
      RETURN
    ENDIF
  ENDIF
  Idp = 0
  RETURN
END SUBROUTINE INDXB
