!** DMACON
SUBROUTINE DMACON
  !>
  !***
  !  Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (MACON-S, DMACON-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  D1MACH
  !***
  ! COMMON BLOCKS    DML5MC

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  USE DML, ONLY : URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  USE service, ONLY : D1MACH
  INTEGER ke
  REAL(8) :: dd
  !* FIRST EXECUTABLE STATEMENT  DMACON
  URO = D1MACH(4)
  SRU = SQRT(URO)
  dd = -LOG10(URO)
  LPAr = INT( 0.5D0*dd )
  ke = INT( 0.5D0 + 0.75D0*dd )
  EPS = 10.0D0**(-2*ke)
  SQOvfl = SQRT(D1MACH(2))
  TWOu = 2.0D0*URO
  FOUru = 4.0D0*URO
END SUBROUTINE DMACON
