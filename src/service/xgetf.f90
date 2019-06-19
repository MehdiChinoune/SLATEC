!** XGETF
SUBROUTINE XGETF(Kontrl)
  !> Return the current value of the error control flag.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XGETF-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !   Abstract
  !        XGETF returns the current value of the error control flag
  !        in KONTRL.  See subroutine XSETF for flag value meanings.
  !        (KONTRL is an output parameter only.)
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

  INTEGER :: Kontrl
  !* FIRST EXECUTABLE STATEMENT  XGETF
  Kontrl = J4SAVE(2,0,.FALSE.)
END SUBROUTINE XGETF
