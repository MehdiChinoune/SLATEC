!** XERDMP
SUBROUTINE XERDMP
  !>
  !  Print the error tables and then clear them.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XERDMP-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XERDMP prints the error tables, then clears them.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  XERSVE

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Changed call of XERSAV to XERSVE.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER kount
  !* FIRST EXECUTABLE STATEMENT  XERDMP
  CALL XERSVE(' ',' ',' ',0,0,0,kount)
END SUBROUTINE XERDMP
