!*==QCDRF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QCDRF
      SUBROUTINE QCDRF(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--QCDRF5
!***BEGIN PROLOGUE  QCDRF
!***PURPOSE  Quick check for DRF.
!***LIBRARY   SLATEC
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Pexton, R. L., (LLNL)
!***DESCRIPTION
!
!            QUICK TEST FOR CARLSON INTEGRAL DRF
!
!***ROUTINES CALLED  D1MACH, DRF, NUMXER, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   910708  Minor modifications in use of KPRINT.  (WRB)
!   930214  Added more digits to ALEM.  (WRB)
!***END PROLOGUE  QCDRF
      INTEGER Kprint , Ipass , contrl , kontrl , Lun , ier
      INTEGER ipass1 , ipass2 , ipass3 , ipass4 , NUMXER
      DOUBLE PRECISION alem , trf , DRF , dif , D1MACH
      EXTERNAL D1MACH , DRF , NUMXER , XERCLR , XGETF , XSETF
!***FIRST EXECUTABLE STATEMENT  QCDRF
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
99001 FORMAT (' DRF - FORCE ERROR 1 TO OCCUR')
      trf = DRF(-1.0D0,-1.0D0,-1.0D0,ier)
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
99002 FORMAT (' DRF - FORCE ERROR 2 TO OCCUR')
      trf = DRF(D1MACH(1),D1MACH(1),D1MACH(1),ier)
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
99003 FORMAT (' DRF - FORCE ERROR 3 TO OCCUR')
      trf = DRF(D1MACH(2),D1MACH(2),D1MACH(2),ier)
      ier = NUMXER(ier)
      IF ( ier==3 ) THEN
        ipass3 = 1
      ELSE
        ipass3 = 0
      ENDIF
      CALL XERCLR
!
!  ARGUMENTS IN RANGE
!  ALEM=LEMNISCATE CONSTANT A
!
      alem = 1.3110287771460599052324197949455597068D0
      trf = DRF(0.0D0,1.0D0,2.0D0,ier)
      CALL XERCLR
      dif = trf - alem
      IF ( (ABS(dif/alem)<1000.0D0*D1MACH(4)).AND.(ier==0) ) THEN
        ipass4 = 1
      ELSE
        ipass4 = 0
      ENDIF
      Ipass = MIN(ipass1,ipass2,ipass3,ipass4)
      IF ( Kprint==0 ) THEN
      ELSEIF ( Kprint==1 ) THEN
        IF ( Ipass/=1 ) WRITE (Lun,99006)
      ELSEIF ( Ipass==1 ) THEN
        WRITE (Lun,99004)
99004   FORMAT (' DRF - PASSED')
      ELSE
        WRITE (Lun,99006)
        IF ( ipass4==0 ) WRITE (Lun,99005) alem , trf , dif
99005   FORMAT (' CORRECT ANSWER =',1PD20.14/'COMPUTED ANSWER =',
     &          D20.14/'     DIFFERENCE =',D20.14)
      ENDIF
      CALL XSETF(contrl)
99006 FORMAT (' DRF - FAILED')
      END SUBROUTINE QCDRF
