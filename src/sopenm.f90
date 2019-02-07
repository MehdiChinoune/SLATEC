!*==SOPENM.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SOPENM
SUBROUTINE SOPENM(Ipage,Lpage)
  IMPLICIT NONE
  !*--SOPENM5
  !*** Start of declarations inserted by SPAG
  INTEGER ios, Ipage, ipagef, Lpage, lpg
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SOPENM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      ALL (SOPENM-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     1. OPEN UNIT NUMBER IPAGEF AS A RANDOM ACCESS FILE.
  !
  !     2. THE RECORD LENGTH IS CONSTANT=LPG.
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !***END PROLOGUE  SOPENM
  CHARACTER(8) :: xern1
  !
  !***FIRST EXECUTABLE STATEMENT  SOPENM
  ipagef = Ipage
  lpg = Lpage
  OPEN (UNIT=ipagef,IOSTAT=ios,ERR=100,STATUS='UNKNOWN',ACCESS='DIRECT',&
    FORM='UNFORMATTED',RECL=lpg)
  RETURN
  !
  100  WRITE (xern1,'(I8)') ios
  CALL XERMSG('SLATEC','SOPENM','IN SPLP, OPEN HAS ERROR FLAG = '//xern1,&
    100,1)
END SUBROUTINE SOPENM
