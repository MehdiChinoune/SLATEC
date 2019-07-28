!** INXCC
PURE SUBROUTINE INXCC(I,Ir,Idxc,Nc)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (INXCC-I)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    CCBLK

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE CCBLK, ONLY : nm_com
  !
  INTEGER, INTENT(IN) :: I, Ir
  INTEGER, INTENT(OUT) :: Idxc, Nc
  !* FIRST EXECUTABLE STATEMENT  INXCC
  Nc = 2**Ir
  Idxc = I
  IF( Idxc+Nc-1>nm_com ) Nc = 0
  !
END SUBROUTINE INXCC