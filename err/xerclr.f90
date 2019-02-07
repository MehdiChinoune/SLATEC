!*==XERCLR.f90  processed by SPAG 6.72Dc at 10:54 on  6 Feb 2019
!DECK XERCLR
SUBROUTINE XERCLR
  IMPLICIT NONE
  !*--XERCLR5
  !*** Start of declarations inserted by SPAG
  INTEGER J4SAVE , junk
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  XERCLR
  !***PURPOSE  Reset current error number to zero.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERCLR-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        This routine simply resets the current error number to zero.
  !        This may be necessary in order to determine that a certain
  !        error has occurred again since the last time NUMXER was
  !        referenced.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  J4SAVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERCLR
  !***FIRST EXECUTABLE STATEMENT  XERCLR
  junk = J4SAVE(1,0,.TRUE.)
END SUBROUTINE XERCLR
