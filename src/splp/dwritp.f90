!DECK DWRITP
SUBROUTINE DWRITP(Ipage,List,Rlist,Lpage,Irec)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DWRITP
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SWRITP-S, DWRITP-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     WRITE RECORD NUMBER IRECN, OF LENGTH LPG, FROM STORAGE
  !     ARRAY LIST(*) ONTO UNIT NUMBER IPAGEF.
  !     WRITE RECORD NUMBER IRECN+1, OF LENGTH LPG, ONTO UNIT
  !     NUMBER IPAGEF FROM THE STORAGE ARRAY RLIST(*).
  !
  !     TO CHANGE THIS PROGRAM UNIT TO DOUBLE PRECISION CHANGE
  !     /REAL (12 BLANKS)/ TO /DOUBLE PRECISION/.
  !
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !***END PROLOGUE  DWRITP
  INTEGER i, Ipage, ipagef, Irec, irecn, Lpage, lpg
  INTEGER List(*)
  REAL(8) :: Rlist(*)
  CHARACTER(8) :: xern1, xern2
  !***FIRST EXECUTABLE STATEMENT  DWRITP
  ipagef = Ipage
  lpg = Lpage
  irecn = Irec
  WRITE (ipagef,REC=irecn,ERR=100) (List(i),i=1,lpg)
  WRITE (ipagef,REC=irecn+1,ERR=100) (Rlist(i),i=1,lpg)
  RETURN
  !
  100  WRITE (xern1,'(I8)') lpg
  WRITE (xern2,'(I8)') irecn
  CALL XERMSG('SLATEC','DWRITP','IN DSPLP, LGP = '//xern1//' IRECN = '//&
    xern2,100,1)
END SUBROUTINE DWRITP