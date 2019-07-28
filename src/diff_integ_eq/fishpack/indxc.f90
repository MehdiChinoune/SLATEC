!** INDXC
PURE SUBROUTINE INDXC(I,Ir,Idxc,Nc)
  !> Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (INDXC-I)
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
  USE CBLKT, ONLY : nm_com
  !
  INTEGER, INTENT(IN) :: I, Ir
  INTEGER, INTENT(OUT) :: Idxc, Nc
  !* FIRST EXECUTABLE STATEMENT  INDXC
  Nc = 2**Ir
  Idxc = I
  IF( Idxc+Nc-1>nm_com ) Nc = 0
  !
END SUBROUTINE INDXC