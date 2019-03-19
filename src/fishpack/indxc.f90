!** INDXC
SUBROUTINE INDXC(I,Ir,Idxc,Nc)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BLKTRI
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
  
  REAL CNV, EPS
  INTEGER I, Idxc, IK, Ir, K, Nc, NCMplx, NM, NPP
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !* FIRST EXECUTABLE STATEMENT  INDXC
  Nc = 2**Ir
  Idxc = I
  IF ( Idxc+Nc-1>NM ) Nc = 0
END SUBROUTINE INDXC
