!** SCLOSM
SUBROUTINE SCLOSM(Ipage)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (SCLOSM-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     1. UNLOAD, RELEASE, OR CLOSE UNIT NUMBER IPAGEF.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  
  INTEGER ios, Ipage, ipagef
  CHARACTER(8) :: xern1
  !
  !* FIRST EXECUTABLE STATEMENT  SCLOSM
  ipagef = Ipage
  CLOSE (UNIT=ipagef,IOSTAT=ios,ERR=100,STATUS='KEEP')
  RETURN
  !
  100  WRITE (xern1,'(I8)') ios
  CALL XERMSG('SLATEC','SCLOSM','IN SPLP, CLOSE HAS ERROR FLAG = '//xern1,100,1)
END SUBROUTINE SCLOSM
