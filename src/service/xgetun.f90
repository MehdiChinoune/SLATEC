!** XGETUN
SUBROUTINE XGETUN(Iunit)
  !>
  !  Return the (first) output file to which error messages
  !            are being sent.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XGETUN-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XGETUN gets the (first) output file to which error messages
  !        are being sent.  To find out if more than one file is being
  !        used, one must use the XGETUA routine.
  !
  !     Description of Parameter
  !      --Output--
  !        IUNIT - the logical unit number of the  (first) unit to
  !                which error messages are being sent.
  !                A value of zero means that the default file, as
  !                defined by the I1MACH routine, is being used.
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

  INTEGER Iunit
  !* FIRST EXECUTABLE STATEMENT  XGETUN
  Iunit = J4SAVE(3,0,.FALSE.)
END SUBROUTINE XGETUN
