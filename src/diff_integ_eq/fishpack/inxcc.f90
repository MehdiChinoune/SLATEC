!** INXCC
SUBROUTINE INXCC(I,Ir,Idxc,Nc)
  USE CCBLK
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBLKTR
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

  INTEGER I, Idxc, Ir, Nc
  !* FIRST EXECUTABLE STATEMENT  INXCC
  Nc = 2**Ir
  Idxc = I
  IF ( Idxc+Nc-1>NM ) Nc = 0
END SUBROUTINE INXCC
