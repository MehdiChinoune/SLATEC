!** XSETUA
SUBROUTINE XSETUA(Iunita,N)
  !>
  !***
  !  Set logical unit numbers (up to 5) to which error
  !            messages are to be sent.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3B
  !***
  ! **Type:**      ALL (XSETUA-A)
  !***
  ! **Keywords:**  ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        XSETUA may be called to declare a list of up to five
  !        logical units, each of which is to receive a copy of
  !        each error message processed by this package.
  !        The purpose of XSETUA is to allow simultaneous printing
  !        of each error message on, say, a main output file,
  !        an interactive terminal, and other files such as graphics
  !        communication files.
  !
  !     Description of Parameters
  !      --Input--
  !        IUNIT - an array of up to five unit numbers.
  !                Normally these numbers should all be different
  !                (but duplicates are not prohibited.)
  !        N     - the number of unit numbers provided in IUNIT
  !                must have 1 .LE. N .LE. 5.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900510  Change call to XERRWV to XERMSG.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, indexx, Iunita(5), junk, N
  CHARACTER(8) :: xern1
  !* FIRST EXECUTABLE STATEMENT  XSETUA
  !
  IF ( N<1.OR.N>5 ) THEN
    WRITE (xern1,'(I8)') N
    CALL XERMSG('SLATEC','XSETUA','INVALID NUMBER OF UNITS, N = '//xern1,1,2)
    RETURN
  END IF
  !
  DO i = 1, N
    indexx = i + 4
    IF ( i==1 ) indexx = 3
    junk = J4SAVE(indexx,Iunita(i),.TRUE.)
  END DO
  junk = J4SAVE(5,N,.TRUE.)
END SUBROUTINE XSETUA
