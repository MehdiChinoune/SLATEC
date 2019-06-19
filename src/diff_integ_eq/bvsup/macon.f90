!** MACON
SUBROUTINE MACON
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (MACON-S, DMACON-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !    Sets up machine constants using R1MACH
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  R1MACH
  !***
  ! COMMON BLOCKS    ML5MCO

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  USE ML, ONLY : uro_com, sru_com, eps_com, sqovfl_com, twou_com, fouru_com, lpar_com
  USE service, ONLY : R1MACH
  REAL(SP) :: dd
  INTEGER :: ke
  !* FIRST EXECUTABLE STATEMENT  MACON
  uro_com = R1MACH(4)
  sru_com = SQRT(uro_com)
  dd = -LOG10(uro_com)
  lpar_com = INT( 0.5*dd )
  ke = INT( 0.5 + 0.75*dd )
  eps_com = 10.**(-2*ke)
  sqovfl_com = SQRT(R1MACH(2))
  twou_com = 2.0*uro_com
  fouru_com = 4.0*uro_com
END SUBROUTINE MACON
