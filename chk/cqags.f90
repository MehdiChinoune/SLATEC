!*==CQAGS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQAGS
SUBROUTINE CQAGS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQAGS5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv, Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQAGS
  !***PURPOSE  Quick check for QAGS.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQAGS-S, CDQAGS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, F0S, F1S, F2S, F3S, F4S, F5S, QAGS, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   911114  Modified test on IER=4 to allow IER=5.  (WRB)
  !***END PROLOGUE  CQAGS
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, error, exact0, &
    exact1, exact2, exact3, exact4, F0S, F1S, F2S, F3S, F4S, &
    F5S, oflow, result, uflow, work
  INTEGER ier, ip, Ipass, iwork, Kprint, last, lenw, limit, neval
  DIMENSION ierv(5), iwork(200), work(800)
  EXTERNAL F0S, F1S, F2S, F3S, F4S, F5S
  DATA exact0/0.2E+01/
  DATA exact1/0.115470066904E+01/
  DATA exact2/0.909864525656E-02/
  DATA exact3/0.31415926535897932E+01/
  DATA exact4/0.19984914554328673E+04/
  !***FIRST EXECUTABLE STATEMENT  CQAGS
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QAGS QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  limit = 200
  lenw = limit*4
  epsabs = 0.0E+00
  epmach = R1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1E-07)
  a = 0.0E+00
  b = 0.1E+01
  CALL QAGS(F0S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  error = ABS(result-exact0)
  ierv(1) = ier
  ip = 0
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL QAGS(F1S,a,b,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,&
    work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 4 OR 1
  !
  uflow = R1MACH(1)
  a = 0.1E+00
  CALL QAGS(F2S,a,b,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 3 OR 4 OR 1 OR 2
  !
  a = 0.0E+00
  b = 0.5E+01
  CALL QAGS(F3S,a,b,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 4, OR 5 OR 3 OR 1 OR 0
  !
  b = 0.1E+01
  epsrel = 1.E-4
  CALL QAGS(F4S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  !      IER=4
  ierv(1) = ier
  ierv(2) = 5
  ierv(3) = 3
  ierv(4) = 1
  ierv(5) = 0
  ip = 0
  IF ( ier==5.OR.ier==4.OR.ier==3.OR.ier==1.OR.ier==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,5)
  !
  ! TEST ON IER = 5
  !
  oflow = R1MACH(2)
  CALL QAGS(F5S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  CALL QAGS(F1S,a,b,epsabs,0.0E+00,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQAGS FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQAGS PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQAGS
