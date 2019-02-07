!*==XERDMP.f90  processed by SPAG 6.72Dc at 10:54 on  6 Feb 2019
!DECK XERDMP
SUBROUTINE XERDMP
  IMPLICIT NONE
  !*--XERDMP5
  !*** Start of declarations inserted by SPAG
  INTEGER kount
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  XERDMP
  !***PURPOSE  Print the error tables and then clear them.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      ALL (XERDMP-A)
  !***KEYWORDS  ERROR, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        XERDMP prints the error tables, then clears them.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  XERSVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Changed call of XERSAV to XERSVE.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  XERDMP
  !***FIRST EXECUTABLE STATEMENT  XERDMP
  CALL XERSVE(' ',' ',' ',0,0,0,kount)
END SUBROUTINE XERDMP
