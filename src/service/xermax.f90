!** XERMAX
SUBROUTINE XERMAX(Maxx)
  IMPLICIT NONE
  !>
  !***
  !  Set maximum number of times any error message is to be
  !            printed.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XERMAX-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XERMAX sets the maximum number of times any message
  !        is to be printed.  That is, non-fatal messages are
  !        not to be printed after they have occurred MAX times.
  !        Such non-fatal messages may be printed less than
  !        MAX times even if they occur MAX times, if error
  !        suppression mode (KONTRL=0) is ever in effect.
  !
  !     Description of Parameter
  !      --Input--
  !        MAX - the maximum number of times any one message
  !              is to be printed.
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

  INTEGER J4SAVE, junk, Maxx
  !* FIRST EXECUTABLE STATEMENT  XERMAX
  junk = J4SAVE(4,Maxx,.TRUE.)
END SUBROUTINE XERMAX
