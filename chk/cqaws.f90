!*==CQAWS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQAWS
SUBROUTINE CQAWS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQAWS5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQAWS
  !***PURPOSE  Quick check for QAWS.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQAWS-S, CDQAWS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, F0WS, F1WS, QAWS, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CQAWS
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  REAL a , abserr , b , R1MACH , epmach , epsabs , epsrel , error , exact0 , &
    exact1 , F0WS , F1WS , alfa , beta , result , uflow , work
  INTEGER ier , ip , Ipass , iwork , Kprint , last , lenw , limit , neval , &
    integr
  DIMENSION work(800) , iwork(200) , ierv(2)
  EXTERNAL F0WS , F1WS
  DATA exact0/0.5350190569223644E+00/
  DATA exact1/0.1998491554328673E+04/
  !***FIRST EXECUTABLE STATEMENT  CQAWS
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWS QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  alfa = -0.5E+00
  beta = -0.5E+00
  integr = 1
  a = 0.0E+00
  b = 0.1E+01
  limit = 200
  lenw = limit*4
  epsabs = 0.0E+00
  epmach = R1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1E-07)
  CALL QAWS(F0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
    limit,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL QAWS(F0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
    2,8,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 1
  !
  uflow = R1MACH(1)
  CALL QAWS(F0WS,a,b,alfa,beta,integr,uflow,0.0E+00,result,abserr,neval,ier,&
    limit,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==2.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 3 OR 1
  !
  CALL QAWS(F1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
    limit,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==3.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 6
  !
  integr = 0
  CALL QAWS(F1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,ier,&
    limit,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQAWS FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQAWS PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQAWS
