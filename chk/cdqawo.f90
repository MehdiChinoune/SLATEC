!*==CDQAWO.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAWO
SUBROUTINE CDQAWO(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CDQAWO5
  !*** Start of declarations inserted by SPAG
  INTEGER leniw
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CDQAWO
  !***PURPOSE  Quick check for DQAWO.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAWO-S, CDQAWO-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF0O, DF1O, DF2O, DPRIN, DQAWO
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CDQAWO
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  DOUBLE PRECISION a , abserr , b , epmach , epsabs , epsrel , error , &
    exact0 , DF0O , DF1O , DF2O , oflow , omega , pi , &
    result , D1MACH , uflow , work
  INTEGER ier , ierv , integr , ip , Ipass , iwork , Kprint , last , lenw , &
    Lun , maxp1 , neval
  DIMENSION work(1325) , iwork(400) , ierv(4)
  EXTERNAL DF0O , DF1O , DF2O
  DATA exact0/0.1042872789432789D+05/
  DATA pi/0.31415926535897932D+01/
  !***FIRST EXECUTABLE STATEMENT  CDQAWO
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWO QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  maxp1 = 21
  leniw = 400
  lenw = leniw*2 + maxp1*25
  epsabs = 0.0D+00
  epmach = D1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1D-07)
  a = 0.0D+00
  b = pi
  omega = 0.1D+01
  integr = 2
  CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  leniw = 2
  lenw = leniw*2 + maxp1*25
  CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 4 OR 1
  !
  uflow = D1MACH(1)
  leniw = 400
  lenw = leniw*2 + maxp1*25
  CALL DQAWO(DF0O,a,b,omega,integr,uflow,0.0D+00,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 3 OR 4 OR 1 OR 2
  !
  b = 0.5D+01
  omega = 0.0D+00
  integr = 1
  CALL DQAWO(DF1O,a,b,omega,integr,uflow,0.0D+00,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 5
  !
  b = 0.1D+01
  oflow = D1MACH(2)
  CALL DQAWO(DF2O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  integr = 3
  CALL DQAWO(DF0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWO FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWO PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAWO
