!DECK CQAWF
SUBROUTINE CQAWF(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER ierv, integr, iwork, leniw, Lun, maxp1
  !***BEGIN PROLOGUE  CQAWF
  !***PURPOSE  Quick check for QAWF.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQAWF-S, CDQAWF-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, F0F, F1F, QAWF, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CQAWF
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL a, abserr, R1MACH, epsabs, epmach, error, exact0, F0F, F1F, &
    omega, pi, result, uflow, work
  INTEGER ier, ip, Ipass, Kprint, lenw, limit, limlst, lst, neval
  DIMENSION ierv(3), iwork(450), work(1425)
  EXTERNAL F0F, F1F
  DATA exact0/0.1422552162575912E+01/
  DATA pi/0.31415926535897932E+01/
  !***FIRST EXECUTABLE STATEMENT  CQAWF
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWF QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  maxp1 = 21
  limlst = 50
  limit = 200
  leniw = limit*2 + limlst
  lenw = leniw*2 + maxp1*25
  epmach = R1MACH(4)
  epsabs = MAX(SQRT(epmach),0.1E-02)
  a = 0.0E+00
  omega = 0.8E+01
  integr = 2
  CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsabs ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  limlst = 3
  leniw = 403
  lenw = leniw*2 + maxp1*25
  CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 3 OR 4 OR 1
  !
  limlst = 50
  leniw = limit*2 + limlst
  lenw = leniw*2 + maxp1*25
  uflow = R1MACH(1)
  CALL QAWF(F1F,a,0.0E+00,1,uflow,result,abserr,neval,ier,limlst,lst,leniw,&
    maxp1,lenw,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 6
  !
  limlst = 50
  leniw = 20
  CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 7
  !
  limlst = 50
  leniw = 52
  lenw = leniw*2 + maxp1*25
  CALL QAWF(F0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==7 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,7,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQAWF FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQAWF PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQAWF
