!DECK NUMXER
FUNCTION NUMXER(Nerr)
  IMPLICIT NONE
  INTEGER J4SAVE, Nerr, NUMXER
  !***BEGIN PROLOGUE  NUMXER
  !***PURPOSE  Return the most recent error number.
  !***LIBRARY   SLATEC (XERROR)
  !***CATEGORY  R3C
  !***TYPE      INTEGER (NUMXER-I)
  !***KEYWORDS  ERROR NUMBER, XERROR
  !***AUTHOR  Jones, R. E., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract
  !        NUMXER returns the most recent error number,
  !        in both NUMXER and the parameter NERR.
  !
  !***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***ROUTINES CALLED  J4SAVE
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910411  Made user-callable and added KEYWORDS section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  NUMXER
  !***FIRST EXECUTABLE STATEMENT  NUMXER
  Nerr = J4SAVE(1,0,.FALSE.)
  NUMXER = Nerr
END FUNCTION NUMXER
