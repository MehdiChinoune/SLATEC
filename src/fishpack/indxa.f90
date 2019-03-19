!** INDXA
SUBROUTINE INDXA(I,Ir,Idxa,Na)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      INTEGER (INDXA-I)
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
  INTEGER I, Idxa, IK, Ir, K, Na, NCMplx, NM, NPP
  COMMON /CBLKT / NPP, K, EPS, CNV, NM, NCMplx, IK
  !* FIRST EXECUTABLE STATEMENT  INDXA
  Na = 2**Ir
  Idxa = I - Na + 1
  IF ( I>NM ) Na = 0
END SUBROUTINE INDXA
