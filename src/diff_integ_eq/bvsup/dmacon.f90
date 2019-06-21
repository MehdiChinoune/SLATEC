!** DMACON
SUBROUTINE DMACON
  !> Subsidiary to DBVSUP
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
  USE DML, ONLY : uro_com, sru_com, eps_com, sqovfl_com, twou_com, fouru_com, lpar_com
  USE service, ONLY : D1MACH
  INTEGER :: ke
  REAL(DP) :: dd
  !* FIRST EXECUTABLE STATEMENT  DMACON
  uro_com = D1MACH(4)
  sru_com = SQRT(uro_com)
  dd = -LOG10(uro_com)
  lpar_com = INT( 0.5_DP*dd )
  ke = INT( 0.5_DP + 0.75_DP*dd )
  eps_com = 10._DP**(-2*ke)
  sqovfl_com = SQRT(D1MACH(2))
  twou_com = 2._DP*uro_com
  fouru_com = 4._DP*uro_com
END SUBROUTINE DMACON
