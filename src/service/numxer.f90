!** NUMXER
INTEGER FUNCTION NUMXER(Nerr)
  IMPLICIT NONE
  !>
  !***
  !  Return the most recent error number.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      INTEGER (NUMXER-I)
  !***
  ! **Keywords:**  ERROR NUMBER, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        NUMXER returns the most recent error number,
  !        in both NUMXER and the parameter NERR.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  J4SAVE

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910411  Made user-callable and added KEYWORDS section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER J4SAVE, Nerr
  !* FIRST EXECUTABLE STATEMENT  NUMXER
  Nerr = J4SAVE(1,0,.FALSE.)
  NUMXER = Nerr
END FUNCTION NUMXER
