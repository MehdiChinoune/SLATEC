!** XSETF
SUBROUTINE XSETF(Kontrl)
  !> Set the error control flag.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3A
  !***
  ! **Type:**      ALL (XSETF-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XSETF sets the error control flag value to KONTRL.
  !        (KONTRL is an input parameter only.)
  !        The following table shows how each message is treated,
  !        depending on the values of KONTRL and LEVEL.  (See XERMSG
  !        for description of LEVEL.)
  !
  !        If KONTRL is zero or negative, no information other than the
  !        message itself (including numeric values, if any) will be
  !        printed.  If KONTRL is positive, introductory messages,
  !        trace-backs, etc., will be printed in addition to the message.
  !
  !              ABS(KONTRL)
  !        LEVEL        0              1              2
  !        value
  !          2        fatal          fatal          fatal
  !
  !          1     not printed      printed         fatal
  !
  !          0     not printed      printed        printed
  !
  !         -1     not printed      printed        printed
  !                                  only           only
  !                                  once           once
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Change call to XERRWV to XERMSG.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: junk, Kontrl
  CHARACTER(8) :: xern1
  !* FIRST EXECUTABLE STATEMENT  XSETF
  IF( ABS(Kontrl)>2 ) THEN
    WRITE (xern1,'(I8)') Kontrl
    CALL XERMSG('XSETF','INVALID ARGUMENT = '//xern1,1,2)
    RETURN
  END IF
  !
  junk = J4SAVE(2,Kontrl,.TRUE.)
END SUBROUTINE XSETF
