!*==CDQAG.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAG
      SUBROUTINE CDQAG(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--CDQAG5
!*** Start of declarations inserted by SPAG
      INTEGER ierv , Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CDQAG
!***PURPOSE  Quick check for DQAG.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (CQAG-S, CDQAG-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  D1MACH, DF1G, DF2G, DF3G, DPRIN, DQAG
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901205  Added PASS/FAIL message and changed the name of the first
!           argument.  (RWC)
!   910501  Added PURPOSE and TYPE records.  (WRB)
!***END PROLOGUE  CDQAG
!
! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
!
      DOUBLE PRECISION a , abserr , b , D1MACH , epmach , epsabs , epsrel , 
     &                 error , exact1 , exact2 , exact3 , DF1G , DF2G , DF3G , 
     &                 pi , result , uflow , work
      INTEGER ier , ip , Ipass , iwork , key , Kprint , last , lenw , limit , 
     &        neval
      DIMENSION ierv(2) , iwork(100) , work(400)
      EXTERNAL DF1G , DF2G , DF3G
      DATA pi/0.31415926535897932D+01/
      DATA exact1/0.1154700538379252D+01/
      DATA exact2/0.11780972450996172D+00/
      DATA exact3/0.1855802D+02/
!***FIRST EXECUTABLE STATEMENT  CDQAG
      IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAG QUICK CHECK''/)')
!
! TEST ON IER = 0
!
      Ipass = 1
      limit = 100
      lenw = limit*4
      epsabs = 0.0D+00
      epmach = D1MACH(4)
      key = 6
      epsrel = MAX(SQRT(epmach),0.1D-07)
      a = 0.0D+00
      b = 0.1D+01
      CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,
     &          last,iwork,work)
      ierv(1) = ier
      ip = 0
      error = ABS(exact1-result)
      IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
!
! TEST ON IER = 1
!
      limit = 1
      lenw = limit*4
      b = pi*0.2D+01
      CALL DQAG(DF2G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,
     &          last,iwork,work)
      ierv(1) = ier
      ip = 0
      IF ( ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,ierv,1)
!
! TEST ON IER = 2 OR 1
!
      uflow = D1MACH(1)
      limit = 100
      lenw = limit*4
      CALL DQAG(DF2G,a,b,uflow,0.0D+00,key,result,abserr,neval,ier,limit,lenw,
     &          last,iwork,work)
      ierv(1) = ier
      ierv(2) = 1
      ip = 0
      IF ( ier==2.OR.ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,2)
!
! TEST ON IER = 3 OR 1
!
      b = 0.1D+01
      CALL DQAG(DF3G,a,b,epsabs,epsrel,1,result,abserr,neval,ier,limit,lenw,
     &          last,iwork,work)
      ierv(1) = ier
      ierv(2) = 1
      ip = 0
      IF ( ier==3.OR.ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,2)
!
! TEST ON IER = 6
!
      lenw = 1
      CALL DQAG(DF1G,a,b,epsabs,epsrel,key,result,abserr,neval,ier,limit,lenw,
     &          last,iwork,work)
      ierv(1) = ier
      ip = 0
      IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.
     &     last==0 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
!
      IF ( Kprint>=1 ) THEN
        IF ( Ipass==0 ) THEN
          WRITE (Lun,'(/'' SOME TEST(S) IN CDQAG FAILED''/)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,'(/'' ALL TEST(S) IN CDQAG PASSED''/)')
        ENDIF
      ENDIF
      END SUBROUTINE CDQAG
