!** XERCLR
SUBROUTINE XERCLR
  !>
  !***
  !  Reset current error number to zero.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XERCLR-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        This routine simply resets the current error number to zero.
  !        This may be necessary in order to determine that a certain
  !        error has occurred again since the last time NUMXER was
  !        referenced.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  J4SAVE

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER junk
  !* FIRST EXECUTABLE STATEMENT  XERCLR
  junk = J4SAVE(1,0,.TRUE.)
END SUBROUTINE XERCLR
