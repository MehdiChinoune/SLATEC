!*==PIMACH.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PIMACH
FUNCTION PIMACH(Dum)
  IMPLICIT NONE
  !*--PIMACH5
  !*** Start of declarations inserted by SPAG
  REAL Dum, PIMACH
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PIMACH
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PIMACH-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This subprogram supplies the value of the constant PI correct to
  !     machine precision where
  !
  !     PI=3.1415926535897932384626433832795028841971693993751058209749446
  !
  !***SEE ALSO  HSTCSP, HSTSSP, HWSCSP
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PIMACH
  !
  !***FIRST EXECUTABLE STATEMENT  PIMACH
  PIMACH = 3.14159265358979
END FUNCTION PIMACH
