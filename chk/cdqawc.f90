!*==CDQAWC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAWC
SUBROUTINE CDQAWC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CDQAWC5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv , Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CDQAWC
  !***PURPOSE  Quick check for DQAWC.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAWC-S, CDQAWC-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF0C, DF1C, DPRIN, DQAWC
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CDQAWC
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL(8) :: a , abserr , b , D1MACH , epmach , epsabs , epsrel , &
    error , exact0 , exact1 , DF0C , DF1C , c , result , &
    uflow , work
  INTEGER ier , ip , Ipass , iwork , Kprint , last , lenw , limit , neval
  DIMENSION work(800) , iwork(200) , ierv(2)
  EXTERNAL DF0C , DF1C
  DATA exact0/ - 0.6284617285065624D+03/
  DATA exact1/0.1855802D+01/
  !***FIRST EXECUTABLE STATEMENT  CDQAWC
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAWC QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  c = 0.5D+00
  a = -1.0D+00
  b = 1.0D+00
  limit = 200
  lenw = limit*4
  epsabs = 0.0D+00
  epmach = D1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1D-07)
  CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ip = 0
  error = ABS(exact0-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact0) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,1,4,last,&
    iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2 OR 1
  !
  uflow = D1MACH(1)
  CALL DQAWC(DF0C,a,b,c,uflow,0.0D+00,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==2.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,2,Kprint,ip,exact0,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 3 OR 1
  !
  CALL DQAWC(DF1C,0.0D+00,b,c,uflow,0.0D+00,result,abserr,neval,ier,limit,&
    lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 1
  ip = 0
  IF ( ier==3.OR.ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,exact1,result,abserr,neval,ierv,2)
  !
  ! TEST ON IER = 6
  !
  epsabs = 0.0D+00
  epsrel = 0.0D+00
  CALL DQAWC(DF0C,a,b,c,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,&
    last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact0,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAWC FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAWC PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAWC
