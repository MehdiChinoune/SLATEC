!*==SWRITP.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SWRITP
      SUBROUTINE SWRITP(Ipage,List,Rlist,Lpage,Irec)
      IMPLICIT NONE
!*--SWRITP5
!*** Start of declarations inserted by SPAG
      INTEGER i , Ipage , ipagef , Irec , irecn , Lpage , lpg
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SWRITP
!***SUBSIDIARY
!***PURPOSE  Subsidiary to SPLP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SWRITP-S, DWRITP-D)
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
!***SEE ALSO  SPLP
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890605  Corrected references to XERRWV.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!***END PROLOGUE  SWRITP
      INTEGER List(*)
      REAL Rlist(*)
      CHARACTER(8) :: xern1 , xern2
!***FIRST EXECUTABLE STATEMENT  SWRITP
      ipagef = Ipage
      lpg = Lpage
      irecn = Irec
      WRITE (ipagef,REC=irecn,ERR=100) (List(i),i=1,lpg)
      WRITE (ipagef,REC=irecn+1,ERR=100) (Rlist(i),i=1,lpg)
      RETURN
!
 100  WRITE (xern1,'(I8)') lpg
      WRITE (xern2,'(I8)') irecn
      CALL XERMSG('SLATEC','SWRITP','IN SPLP, LGP = '//xern1//' IRECN = '//
     &            xern2,100,1)
      END SUBROUTINE SWRITP
