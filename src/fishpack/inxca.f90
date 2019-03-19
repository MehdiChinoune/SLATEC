!** INXCA
SUBROUTINE INXCA(I,Ir,Idxa,Na)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (INXCA-I)
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
  INTEGER I, Idxa, IK, Ir, K, Na, NCMplx, NM, NPP
  COMMON /CCBLK / NPP, K, EPS, CNV, NM, NCMplx, IK
  !* FIRST EXECUTABLE STATEMENT  INXCA
  Na = 2**Ir
  Idxa = I - Na + 1
  IF ( I>NM ) Na = 0
END SUBROUTINE INXCA
