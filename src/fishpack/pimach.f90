!** PIMACH
REAL FUNCTION PIMACH(Dum)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to HSTCSP, HSTSSP and HWSCSP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PIMACH-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subprogram supplies the value of the constant PI correct to
  !     machine precision where
  !
  !     PI=3.1415926535897932384626433832795028841971693993751058209749446
  !
  !***
  ! **See also:**  HSTCSP, HSTSSP, HWSCSP
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL Dum
  !* FIRST EXECUTABLE STATEMENT  PIMACH
  PIMACH = 3.14159265358979
END FUNCTION PIMACH
