!*==CDQAWS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAWS
SUBROUTINE CDQAWS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CDQAWS5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv, Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CDQAWS
  !***PURPOSE  Quick check for DQAWS.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAWS-S, CDQAWS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF0WS, DF1WS, DPRIN, DQAWS
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CDQAWS
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
    error, exact0, exact1, DF0WS, DF1WS, alfa, beta, &
    result, uflow, work
  INTEGER ier, ip, Ipass, iwork, Kprint, last, lenw, limit, neval, &
    integr
  DIMENSION work(800), iwork(200), ierv(2)
  EXTERNAL DF0WS, DF1WS
  DATA exact0/0.5350190569223644D+00/
  DATA exact1/0.1998491554328673D+04/
  !***FIRST EXECUTABLE STATEMENT  CDQAWS
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWS QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  alfa = -0.5D+00
  beta = -0.5D+00
  integr = 1
  a = 0.0D+00
  b = 0.1D+01
  limit = 200
  lenw = limit*4
  epsabs = 0.0D+00
  epmach = D1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1D-07)
  CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
    ier,limit,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL DQAWS(DF0WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
    ier,2,8,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 1
  !
  uflow = D1MACH(1)
  CALL DQAWS(DF0WS,a,b,alfa,beta,integr,uflow,0.0D+00,result,abserr,neval,&
    ier,limit,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==2.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 3 OR 1
  !
  CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
    ier,limit,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==3.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 6
  !
  integr = 0
  CALL DQAWS(DF1WS,a,b,alfa,beta,integr,epsabs,epsrel,result,abserr,neval,&
    ier,limit,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWS FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWS PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAWS
