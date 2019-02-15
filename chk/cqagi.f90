!DECK CQAGI
SUBROUTINE CQAGI(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER ierv, inf
  !***BEGIN PROLOGUE  CQAGI
  !***PURPOSE  Quick check for QAGI.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQAGI-S, CDQAGI-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, QAGI, R1MACH, T0, T1, T2, T3, T4, T5
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891009  Removed unreferenced variables.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CQAGI
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL abserr, bound, R1MACH, epmach, epsabs, epsrel, error, exact0, &
    exact1, exact2, exact3, exact4, oflow, result, T0, T1, T2, &
    T3, T4, T5, uflow, work
  INTEGER ier, ip, Ipass, iwork, Kprint, last, lenw, limit, Lun, &
    neval
  DIMENSION work(800), iwork(200), ierv(4)
  EXTERNAL T0, T1, T2, T3, T4, T5
  DATA exact0/2.0E+00/, exact1/0.115470066904E1/
  DATA exact2/0.909864525656E-02/
  DATA exact3/0.31415926535897932E+01/
  DATA exact4/0.19984914554328673E+04/
  !***FIRST EXECUTABLE STATEMENT  CQAGI
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QAGI QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  limit = 200
  lenw = limit*4
  epsabs = 0.0E+00
  epmach = R1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1E-07)
  bound = 0.0E+00
  inf = 1
  CALL QAGI(T0,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  error = ABS(result-exact0)
  ierv(1) = ier
  ip = 0
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL QAGI(T1,bound,inf,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
    iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 4 OR 1
  !
  uflow = R1MACH(1)
  CALL QAGI(T2,bound,inf,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 3 OR 4 OR 1
  !
  CALL QAGI(T3,bound,inf,uflow,0.0E+00,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 4 OR 3 OR 1
  !
  CALL QAGI(T4,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ierv(2) = 3
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==4.OR.ier==3.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 5
  !
  oflow = R1MACH(2)
  CALL QAGI(T5,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  CALL QAGI(T1,bound,inf,epsabs,0.0E+00,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQAGI FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQAGI PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQAGI
