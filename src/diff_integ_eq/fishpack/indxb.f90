!** INDXB
SUBROUTINE INDXB(I,Ir,Idx,Idp)
  !> Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (INDXB-I)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  BLKTRI
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CBLKT

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920422  Added statement so IDX would always be defined.  (WRB)
  USE CBLKT, ONLY : ik_com, nm_com
  INTEGER :: I, Idp, Idx, Ir
  INTEGER :: id, ipl, izh
  !* FIRST EXECUTABLE STATEMENT  INDXB
  Idx = max(I,1)
  Idp = 0
  IF( Ir<0 ) RETURN
  IF( Ir==0 ) THEN
    IF( I>nm_com ) RETURN
    Idx = I
    Idp = 1
    RETURN
  ELSE
    izh = 2**Ir
    id = I - izh - izh
    Idx = id + id + (Ir-1)*ik_com + Ir + (ik_com-I)/izh + 4
    ipl = izh - 1
    Idp = izh + izh - 1
    IF( I-ipl<=nm_com ) THEN
      IF( I+ipl>nm_com ) Idp = nm_com + ipl - I + 1
      RETURN
    END IF
  END IF
  Idp = 0
  RETURN
END SUBROUTINE INDXB
