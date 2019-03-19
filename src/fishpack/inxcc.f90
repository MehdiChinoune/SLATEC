!** INXCC
SUBROUTINE INXCC(I,Ir,Idxc,Nc)
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
  
  REAL CNV, EPS
  INTEGER I, Idxc, IK, Ir, K, Nc, NCMplx, NM, NPP
  COMMON /CCBLK / NPP, K, EPS, CNV, NM, NCMplx, IK
  !* FIRST EXECUTABLE STATEMENT  INXCC
  Nc = 2**Ir
  Idxc = I
  IF ( Idxc+Nc-1>NM ) Nc = 0
END SUBROUTINE INXCC
