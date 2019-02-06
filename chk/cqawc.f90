!*==CQAWC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQAWC
      SUBROUTINE CQAWC(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--CQAWC5
!*** Start of declarations inserted by SPAG
      INTEGER ierv , Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CQAWC
!***PURPOSE  Quick check for QAWC.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CQAWC-S, CDQAWC-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  CPRIN, F0C, F1C, QAWC, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901205  Added PASS/FAIL message and changed the name of the first
!           argument.  (RWC)
!   910501  Added PURPOSE and TYPE records.  (WRB)
!***END PROLOGUE  CQAWC
!
! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
!
      REAL a , abserr , b , R1MACH , epmach , epsabs , epsrel , error , exact0 , 
     &     exact1 , F0C , F1C , c , result , uflow , work
      INTEGER ier , ip , Ipass , iwork , Kprint , last , lenw , limit , neval
      DIMENSION work(800) , iwork(200) , ierv(2)
      EXTERNAL F0C , F1C
      DATA exact0/ - 0.6284617285065624E+03/
      DATA exact1/0.1855802E+01/
!***FIRST EXECUTABLE STATEMENT  CQAWC
      IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWC QUICK CHECK''/)')
!
! TEST ON IER = 0
!
      Ipass = 1
      c = 0.5E+00
      a = -1.0E+00
      b = 1.0E+00
      limit = 200
      lenw = limit*4
      epsabs = 0.0E+00
      epmach = R1MACH(4)
      epsrel = MAX(SQRT(epmach),0.1E-07)
      CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,
     &          iwork,work)
      ierv(1) = ier
      ip = 0
      error = ABS(exact0-result)
      IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
!
! TEST ON IER = 1
!
      CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,
     &          work)
      ierv(1) = ier
      ip = 0
      IF ( ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
!
! TEST ON IER = 2 OR 1
!
      uflow = R1MACH(1)
      CALL QAWC(F0C,a,b,c,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,
     &          iwork,work)
      ierv(1) = ier
      ierv(2) = 1
      ip = 0
      IF ( ier==2.OR.ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
!
! TEST ON IER = 3 OR 1
!
      CALL QAWC(F1C,0.0E+00,b,c,uflow,0.0E+00,result,abserr,neval,ier,limit,
     &          lenw,last,iwork,work)
      ierv(1) = ier
      ierv(2) = 1
      ip = 0
      IF ( ier==3.OR.ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL CPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
!
! TEST ON IER = 6
!
      epsabs = 0.0E+00
      epsrel = 0.0E+00
      CALL QAWC(F0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,
     &          iwork,work)
      ierv(1) = ier
      ip = 0
      IF ( ier==6 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
!
      IF ( Kprint>=1 ) THEN
        IF ( Ipass==0 ) THEN
          WRITE (Lun,'(/'' SOME TEST(S) IN CQAWC FAILED''/)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,'(/'' ALL TEST(S) IN CQAWC PASSED''/)')
        ENDIF
      ENDIF
      END SUBROUTINE CQAWC
