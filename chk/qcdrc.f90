!*==QCDRC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QCDRC
      SUBROUTINE QCDRC(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--QCDRC5
!***BEGIN PROLOGUE  QCDRC
!***PURPOSE  Quick check for DRC.
!***LIBRARY   SLATEC
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Pexton, R. L., (LLNL)
!***DESCRIPTION
!
!            QUICK TEST FOR CARLSON INTEGRAL DRC
!
!***ROUTINES CALLED  D1MACH, DRC, NUMXER, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   910708  Minor modifications in use of KPRINT.  (WRB)
!***END PROLOGUE  QCDRC
      INTEGER Kprint , Ipass , contrl , kontrl , Lun , ier
      INTEGER ipass1 , ipass2 , ipass3 , ipass4 , NUMXER
      DOUBLE PRECISION pi , trc , DRC , dif , D1MACH
      EXTERNAL D1MACH , DRC , NUMXER , XERCLR , XGETF , XSETF
!***FIRST EXECUTABLE STATEMENT  QCDRC
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
99001 FORMAT (' DRC - FORCE ERROR 1 TO OCCUR')
      trc = DRC(-1.0D0,-1.0D0,ier)
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
99002 FORMAT (' DRC - FORCE ERROR 2 TO OCCUR')
      trc = DRC(D1MACH(1),D1MACH(1),ier)
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
99003 FORMAT (' DRC - FORCE ERROR 3 TO OCCUR')
      trc = DRC(D1MACH(2),D1MACH(2),ier)
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
      pi = 3.141592653589793238462643383279D0
      trc = DRC(0.0D0,0.25D0,ier)
      CALL XERCLR
      dif = trc - pi
      IF ( (ABS(dif/pi)<1000.0D0*D1MACH(4)).AND.(ier==0) ) THEN
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
99004   FORMAT (' DRC - PASSED')
      ELSE
        WRITE (Lun,99006)
        IF ( ipass4==0 ) WRITE (Lun,99005) pi , trc , dif
99005   FORMAT (' CORRECT ANSWER =',1PD20.14/'COMPUTED ANSWER =',
     &          D20.14/'     DIFFERENCE =',D20.14)
      ENDIF
      CALL XSETF(contrl)
99006 FORMAT (' DRC - FAILED')
      END SUBROUTINE QCDRC
