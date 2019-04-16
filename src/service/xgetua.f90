!** XGETUA
SUBROUTINE XGETUA(Iunita,N)
  !>
  !***
  !  Return unit number(s) to which error messages are being
  !            sent.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XGETUA-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XGETUA may be called to determine the unit number or numbers
  !        to which error messages are being sent.
  !        These unit numbers may have been set by a call to XSETUN,
  !        or a call to XSETUA, or may be a default value.
  !
  !     Description of Parameters
  !      --Output--
  !        IUNIT - an array of one to five unit numbers, depending
  !                on the value of N.  A value of zero refers to the
  !                default unit, as defined by the I1MACH machine
  !                constant routine.  Only IUNIT(1),...,IUNIT(N) are
  !                defined by XGETUA.  The values of IUNIT(N+1),...,
  !                IUNIT(5) are not defined (for N .LT. 5) or altered
  !                in any way by XGETUA.
  !        N     - the number of units to which copies of the
  !                error messages are being sent.  N will be in the
  !                range from 1 to 5.
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

  INTEGER i, indexx, Iunita(5), N
  !* FIRST EXECUTABLE STATEMENT  XGETUA
  N = J4SAVE(5,0,.FALSE.)
  DO i = 1, N
    indexx = i + 4
    IF ( i==1 ) indexx = 3
    Iunita(i) = J4SAVE(indexx,0,.FALSE.)
  END DO
END SUBROUTINE XGETUA
