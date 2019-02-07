!*==CDQAWF.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAWF
SUBROUTINE CDQAWF(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CDQAWF5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv, integr, iwork, leniw, Lun, maxp1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CDQAWF
  !***PURPOSE  Quick check for DQAWF.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAWF-S, CDQAWF-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF0F, DF1F, DPRIN, DQAWF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CDQAWF
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL(8) :: a, abserr, D1MACH, epsabs, epmach, error, exact0, &
    DF0F, DF1F, omega, pi, result, uflow, work
  INTEGER ier, ip, Ipass, Kprint, lenw, limit, limlst, lst, neval
  DIMENSION ierv(4), iwork(450), work(1425)
  EXTERNAL DF0F, DF1F
  DATA exact0/0.1422552162575912D+01/
  DATA pi/0.31415926535897932D+01/
  !***FIRST EXECUTABLE STATEMENT  CDQAWF
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWF QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  maxp1 = 21
  limlst = 50
  limit = 200
  leniw = limit*2 + limlst
  lenw = leniw*2 + maxp1*25
  epmach = D1MACH(4)
  epsabs = MAX(SQRT(epmach),0.1D-02)
  a = 0.0D+00
  omega = 0.8D+01
  integr = 2
  CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsabs ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  limlst = 3
  leniw = 403
  lenw = leniw*2 + maxp1*25
  CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 3 OR 4 OR 1 OR 2
  !
  limlst = 50
  leniw = limit*2 + limlst
  lenw = leniw*2 + maxp1*25
  uflow = D1MACH(1)
  CALL DQAWF(DF1F,a,0.0D+00,1,uflow,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 6
  !
  limlst = 50
  leniw = 20
  CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 7
  !
  limlst = 50
  leniw = 52
  lenw = leniw*2 + maxp1*25
  CALL DQAWF(DF0F,a,omega,integr,epsabs,result,abserr,neval,ier,limlst,lst,&
    leniw,maxp1,lenw,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==7 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,7,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWF FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWF PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAWF
