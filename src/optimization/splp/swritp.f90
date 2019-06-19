!** SWRITP
SUBROUTINE SWRITP(Ipage,List,Rlist,Lpage,Irec)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SWRITP-S, DWRITP-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     WRITE RECORD NUMBER IRECN, OF LENGTH LPG, FROM STORAGE
  !     ARRAY LIST(*) ONTO UNIT NUMBER IPAGEF.
  !     WRITE RECORD NUMBER IRECN+1, OF LENGTH LPG, ONTO UNIT
  !     NUMBER IPAGEF FROM THE STORAGE ARRAY RLIST(*).
  !
  !     TO CHANGE THIS PROGRAM UNIT TO DOUBLE PRECISION CHANGE
  !     /REAL (12 BLANKS)/ TO /DOUBLE PRECISION/.
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  USE service, ONLY : XERMSG
  INTEGER :: Ipage, Irec, Lpage
  INTEGER :: List(Lpage)
  REAL(SP) :: Rlist(Lpage)
  INTEGER :: i, ipagef, irecn, lpg
  CHARACTER(8) :: xern1, xern2
  !* FIRST EXECUTABLE STATEMENT  SWRITP
  ipagef = Ipage
  lpg = Lpage
  irecn = Irec
  WRITE (ipagef,REC=irecn,ERR=100) (List(i),i=1,lpg)
  WRITE (ipagef,REC=irecn+1,ERR=100) (Rlist(i),i=1,lpg)
  RETURN
  !
  100  WRITE (xern1,'(I8)') lpg
  WRITE (xern2,'(I8)') irecn
  CALL XERMSG('SWRITP','IN SPLP, LGP = '//xern1//' IRECN = '//&
    xern2,100,1)
END SUBROUTINE SWRITP
