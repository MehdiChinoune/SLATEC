!*==CDQAGP.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQAGP
SUBROUTINE CDQAGP(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CDQAGP5
  !*** Start of declarations inserted by SPAG
  INTEGER ierv
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CDQAGP
  !***PURPOSE  Quick check for DQAGP.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (CQAGP-S, CDQAGP-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DF1P, DF2P, DF3P, DF4P, DPRIN, DQAGP
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CDQAGP
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  DOUBLE PRECISION a , abserr , b , D1MACH , epmach , epsabs , epsrel , &
    error , exact1 , exact2 , exact3 , DF1P , DF2P , DF3P , &
    DF4P , oflow , points , p1 , p2 , result , uflow , work
  INTEGER ier , ip , Ipass , iwork , Kprint , last , leniw , lenw , limit , &
    Lun , neval , npts2
  DIMENSION ierv(4) , iwork(205) , points(5) , work(405)
  EXTERNAL DF1P , DF2P , DF3P , DF4P
  DATA exact1/0.4285277667368085D+01/
  DATA exact2/0.909864525656D-2/
  DATA exact3/0.31415926535897932D+01/
  DATA p1/0.1428571428571428D+00/
  DATA p2/0.6666666666666667D+00/
  !***FIRST EXECUTABLE STATEMENT  CDQAGP
  IF ( Kprint>=2 ) WRITE (Lun,'(''1DQAGP QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  npts2 = 4
  limit = 100
  leniw = limit*2 + npts2
  lenw = limit*4 + npts2
  epsabs = 0.0D+00
  epmach = D1MACH(4)
  epsrel = MAX(SQRT(epmach),0.1D-07)
  a = 0.0D+00
  b = 0.1D+01
  points(1) = p1
  points(2) = p2
  CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  error = ABS(result-exact1)
  ierv(1) = ier
  ip = 0
  IF ( ier==0.AND.error<=epsrel*ABS(exact1) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 1
  !
  leniw = 10
  lenw = leniw*2 - npts2
  CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,1,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 2, 4, 1 OR 3
  !
  npts2 = 3
  points(1) = 0.1D+00
  leniw = limit*2 + npts2
  lenw = limit*4 + npts2
  uflow = D1MACH(1)
  a = 0.1D+00
  CALL DQAGP(DF2P,a,b,npts2,points,uflow,0.0D+00,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 3
  ip = 0
  IF ( ier==2.OR.ier==4.OR.ier==1.OR.ier==3 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,2,Kprint,ip,exact2,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 3 OR 4 OR 1 OR 2
  !
  npts2 = 2
  leniw = limit*2 + npts2
  lenw = limit*4 + npts2
  a = 0.0D+00
  b = 0.5D+01
  CALL DQAGP(DF3P,a,b,npts2,points,uflow,0.0D+00,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  ierv(1) = ier
  ierv(2) = 4
  ierv(3) = 1
  ierv(4) = 2
  ip = 0
  IF ( ier==3.OR.ier==4.OR.ier==1.OR.ier==2 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,3,Kprint,ip,exact3,result,abserr,neval,ierv,4)
  !
  ! TEST ON IER = 5
  !
  b = 0.1D+01
  CALL DQAGP(DF4P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==5 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  oflow = D1MACH(2)
  CALL DPRIN(Lun,5,Kprint,ip,oflow,result,abserr,neval,ierv,1)
  !
  ! TEST ON IER = 6
  !
  npts2 = 5
  leniw = limit*2 + npts2
  lenw = limit*4 + npts2
  points(1) = p1
  points(2) = p2
  points(3) = 0.3D+01
  CALL DQAGP(DF1P,a,b,npts2,points,epsabs,epsrel,result,abserr,neval,ier,&
    leniw,lenw,last,iwork,work)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0.AND.&
    last==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CDQAGP FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CDQAGP PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CDQAGP
