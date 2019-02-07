!*==INDXB.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK INDXB
SUBROUTINE INDXB(I,Ir,Idx,Idp)
  IMPLICIT NONE
  !*--INDXB5
  !*** Start of declarations inserted by SPAG
  REAL CNV , EPS
  INTEGER I , id , Idp , Idx , IK , ipl , Ir , izh , K , NCMplx , NM , NPP
  !*** End of declarations inserted by SPAG
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
  COMMON /CBLKT / NPP , K , EPS , CNV , NM , NCMplx , IK
  !***FIRST EXECUTABLE STATEMENT  INDXB
  Idx = I
  Idp = 0
  IF ( Ir<0 ) GOTO 99999
  IF ( Ir==0 ) THEN
    IF ( I>NM ) GOTO 99999
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
      GOTO 99999
    ENDIF
  ENDIF
  Idp = 0
  RETURN
  99999 END SUBROUTINE INDXB
