!*==CQAWO.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQAWO
SUBROUTINE CQAWO(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQAWO5
  !*** Start of declarations inserted by SPAG
  INTEGER leniw
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQAWO
  !***PURPOSE  Quick check for QAWO.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQAWO-S, CDQAWO-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, F0O, F1O, F2O, QAWO, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CQAWO
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL a , abserr , b , epmach , epsabs , epsrel , error , exact0 , F0O , &
    F1O , F2O , oflow , omega , pi , result , R1MACH , uflow , work
  INTEGER ier , ierv , integr , ip , Ipass , iwork , Kprint , last , lenw , &
    Lun , maxp1 , neval
  DIMENSION work(1325) , iwork(400) , ierv(4)
  EXTERNAL F0O , F1O , F2O
  DATA exact0/0.1042872789432789E+05/
  DATA pi/0.31415926535897932E+01/
  !***FIRST EXECUTABLE STATEMENT  CQAWO
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QAWO QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  maxp1 = 21
  leniw = 400
  lenw = leniw*2 + maxp1*25
  epsabs = 0.0E+00
  epmach = R1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1E-07)
  a = 0.0E+00
  b = pi
  omega = 0.1E+01
  integr = 2
  CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  leniw = 2
  lenw = leniw*2 + maxp1*25
  CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 4 OR 1
  !
  uflow = R1MACH(1)
  leniw = 400
  lenw = leniw*2 + maxp1*25
  CALL QAWO(F0O,a,b,omega,integr,uflow,0.0E+00,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 3 OR 4 OR 1
  !
  b = 0.5E+01
  omega = 0.0E+00
  integr = 1
  CALL QAWO(F1O,a,b,omega,integr,uflow,0.0E+00,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,3,Kprint,ip,pi,result,abserr,neval,ierv,3)
  !
  ! TEST ON IER = 5
  !
  b = 0.1E+01
  oflow = R1MACH(2)
  CALL QAWO(F2O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  integr = 3
  CALL QAWO(F0O,a,b,omega,integr,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,maxp1,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL CPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQAWO FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQAWO PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQAWO
