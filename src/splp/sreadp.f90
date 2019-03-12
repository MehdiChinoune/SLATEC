!DECK SREADP
SUBROUTINE SREADP(Ipage,List,Rlist,Lpage,Irec)
  IMPLICIT NONE
  INTEGER i, Ipage, ipagef, Irec, irecn, Lpage, lpg
  !***BEGIN PROLOGUE  SREADP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SPLP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SREADP-S, DREADP-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     READ RECORD NUMBER IRECN, OF LENGTH LPG, FROM UNIT
  !     NUMBER IPAGEF INTO THE STORAGE ARRAY, LIST(*).
  !     READ RECORD  IRECN+1, OF LENGTH LPG, FROM UNIT NUMBER
  !     IPAGEF INTO THE STORAGE ARRAY RLIST(*).
  !
  !     TO CONVERT THIS PROGRAM UNIT TO DOUBLE PRECISION CHANGE
  !     /REAL (12 BLANKS)/ TO /DOUBLE PRECISION/.
  !
  !***SEE ALSO  SPLP
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !***END PROLOGUE  SREADP
  INTEGER List(*)
  REAL Rlist(*)
  CHARACTER(8) :: xern1, xern2
  !***FIRST EXECUTABLE STATEMENT  SREADP
  ipagef = Ipage
  lpg = Lpage
  irecn = Irec
  READ (ipagef,REC=irecn,ERR=100) (List(i),i=1,lpg)
  READ (ipagef,REC=irecn+1,ERR=100) (Rlist(i),i=1,lpg)
  RETURN
  !
  100  WRITE (xern1,'(I8)') lpg
  WRITE (xern2,'(I8)') irecn
  CALL XERMSG('SLATEC','SREADP','IN SPLP, LPG = '//xern1//' IRECN = '//&
    xern2,100,1)
END SUBROUTINE SREADP
