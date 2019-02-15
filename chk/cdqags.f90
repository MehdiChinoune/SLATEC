!DECK CDQAGS
SUBROUTINE CDQAGS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER ierv, Lun
  !***BEGIN PROLOGUE  CDQAGS
  !***PURPOSE  Quick check for DQAGS.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAGS-S, CDQAGS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF0S, DF1S, DF2S, DF3S, DF4S, DF5S, DPRIN,
  !                    DQAGS
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   911114  Modified test on IER=4 to allow IER=5.  (WRB)
  !***END PROLOGUE  CDQAGS
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL(8) :: a, abserr, b, D1MACH, epmach, epsabs, epsrel, &
    error, exact0, exact1, exact2, exact3, exact4, &
    DF0S, DF1S, DF2S, DF3S, DF4S, DF5S, oflow, &
    result, uflow, work
  INTEGER ier, ip, Ipass, iwork, Kprint, last, lenw, limit, neval
  DIMENSION ierv(5), iwork(200), work(800)
  EXTERNAL DF0S, DF1S, DF2S, DF3S, DF4S, DF5S
  DATA exact0/0.2D+01/
  DATA exact1/0.115470066904D+01/
  DATA exact2/0.909864525656D-02/
  DATA exact3/0.31415926535897932D+01/
  DATA exact4/0.19984914554328673D+04/
  !***FIRST EXECUTABLE STATEMENT  CDQAGS
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAGS QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  limit = 200
  lenw = limit*4
  epsabs = 0.0D+00
  epmach = D1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1D-07)
  a = 0.0D+00
  b = 0.1D+01
  CALL DQAGS(DF0S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  error = ABS(result-exact0)
  ierv(1) = ier
  ip = 0
  IF ( ier==0.AND.error<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL DQAGS(DF1S,a,b,epsabs,epsrel,result,abserr,neval,ier,1,4,last,iwork,&
    work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 4 OR 1
  !
  uflow = D1MACH(1)
  a = 0.1D+00
  CALL DQAGS(DF2S,a,b,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 3 OR 4 OR 1 OR 2
  !
  a = 0.0D+00
  b = 0.5D+01
  CALL DQAGS(DF3S,a,b,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 4, OR 5 OR 3 OR 1 OR 0
  !
  b = 0.1D+01
  CALL DQAGS(DF4S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ierv(2) = 5
  ierv(3) = 3
  ierv(4) = 1
  ierv(5) = 0
  ip = 0
  IF ( ier==5.OR.ier==4.OR.ier==3.OR.ier==1.OR.ier==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,4,Kprint,ip,exact4,result,abserr,neval,ierv,5)
  !
  ! TEST ON IER = 5
  !
  oflow = D1MACH(2)
  CALL DQAGS(DF5S,a,b,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,&
    iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  CALL DQAGS(DF1S,a,b,epsabs,0.0D+00,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGS FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGS PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAGS
