!** XERHLT
SUBROUTINE XERHLT(Messg)
  !>
  !***
  !  Abort program execution and print error message.
  !***
  ! **Library:**   SLATEC (XERROR)
  !***
  ! **Category:**  R3C
  !***
  ! **Type:**      ALL (XERHLT-A)
  !***
  ! **Keywords:**  ABORT PROGRAM EXECUTION, ERROR, XERROR
  !***
  ! **Author:**  Jones, R. E., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !        ***Note*** machine dependent routine
  !        XERHLT aborts the execution of the program.
  !        The error message causing the abort is given in the calling
  !        sequence, in case one needs it for printing on a dayfile,
  !        for example.
  !
  !     Description of Parameters
  !        MESSG is as in XERMSG.
  !
  !***
  ! **References:**  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
  !                 Error-handling Package, SAND82-0800, Sandia
  !                 Laboratories, 1982.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900206  Routine changed from user-callable to subsidiary.  (WRB)
  !   900510  Changed calling sequence to delete length of character
  !           and changed routine name from XERABT to XERHLT.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  CHARACTER*(*) Messg
  !* FIRST EXECUTABLE STATEMENT  XERHLT
  STOP
END SUBROUTINE XERHLT
