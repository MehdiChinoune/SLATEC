!DECK QCRC
SUBROUTINE QCRC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QCRC
  !***PURPOSE  Quick check for RC.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Pexton, R. L., (LLNL)
  !***DESCRIPTION
  !
  !            QUICK TEST FOR CARLSON INTEGRAL RC
  !
  !***ROUTINES CALLED  NUMXER, R1MACH, RC, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   790801  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !***END PROLOGUE  QCRC
  INTEGER Kprint, Ipass, contrl, kontrl, Lun, ier
  INTEGER ipass1, ipass2, ipass3, ipass4, NUMXER
  REAL pi, trc, RC, dif, R1MACH
  EXTERNAL NUMXER, R1MACH, RC, XERCLR, XGETF, XSETF
  !***FIRST EXECUTABLE STATEMENT  QCRC
  CALL XERCLR
  CALL XGETF(contrl)
  IF ( Kprint>=3 ) THEN
    kontrl = +1
  ELSE
    kontrl = 0
  ENDIF
  CALL XSETF(kontrl)
  !
  !  FORCE ERROR 1
  !
  IF ( Kprint>=3 ) WRITE (Lun,99001)
  99001 FORMAT (' RC - FORCE ERROR 1 TO OCCUR')
  trc = RC(-1.0E0,-1.0E0,ier)
  ier = NUMXER(ier)
  IF ( ier==1 ) THEN
    ipass1 = 1
  ELSE
    ipass1 = 0
  ENDIF
  CALL XERCLR
  !
  !  FORCE ERROR 2
  !
  IF ( Kprint>=3 ) WRITE (Lun,99002)
  99002 FORMAT (' RC - FORCE ERROR 2 TO OCCUR')
  trc = RC(R1MACH(1),R1MACH(1),ier)
  ier = NUMXER(ier)
  IF ( ier==2 ) THEN
    ipass2 = 1
  ELSE
    ipass2 = 0
  ENDIF
  CALL XERCLR
  !
  !  FORCE ERROR 3
  !
  IF ( Kprint>=3 ) WRITE (Lun,99003)
  99003 FORMAT (' RC - FORCE ERROR 3 TO OCCUR')
  trc = RC(R1MACH(2),R1MACH(2),ier)
  ier = NUMXER(ier)
  IF ( ier==3 ) THEN
    ipass3 = 1
  ELSE
    ipass3 = 0
  ENDIF
  CALL XERCLR
  !
  !  ARGUMENTS IN RANGE
  !
  pi = 3.1415926535897932E0
  trc = RC(0.0E0,0.25E0,ier)
  CALL XERCLR
  dif = trc - pi
  IF ( (ABS(dif/pi)<1000.0E0*R1MACH(4)).AND.(ier==0) ) THEN
    ipass4 = 1
  ELSE
    ipass4 = 0
  ENDIF
  Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
  IF ( Kprint<=0 ) THEN
  ELSEIF ( Kprint==1 ) THEN
    IF ( Ipass/=1 ) WRITE (Lun,99006)
  ELSEIF ( Ipass==1 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (' RC - PASSED')
  ELSE
    WRITE (Lun,99006)
    IF ( ipass4==0 ) WRITE (Lun,99005) pi, trc, dif
    99005   FORMAT (' CORRECT ANSWER =',1PE14.6/'COMPUTED ANSWER =',&
      E14.6/'     DIFFERENCE =',E14.6)
  ENDIF
  CALL XSETF(contrl)
  99006 FORMAT (' RC - FAILED')
END SUBROUTINE QCRC
