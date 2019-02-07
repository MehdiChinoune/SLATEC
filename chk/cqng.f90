!*==CQNG.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CQNG
SUBROUTINE CQNG(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--CQNG5
  !*** Start of declarations inserted by SPAG
  INTEGER Lun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  CQNG
  !***PURPOSE  Quick check for QNG.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CQNG-S, CDQNG-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  CPRIN, F1N, F2N, QNG, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Added PASS/FAIL message and changed the name of the first
  !           argument.  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !***END PROLOGUE  CQNG
  !
  ! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
  !
  REAL a, abserr, b, R1MACH, epmach, epsabs, epsrel, exact1, error, &
    exact2, F1N, F2N, result, uflow
  INTEGER ier, ierv, ip, Ipass, Kprint, neval
  DIMENSION ierv(1)
  EXTERNAL F1N, F2N
  DATA exact1/0.7281029132255818E+00/
  DATA exact2/0.1E+02/
  !***FIRST EXECUTABLE STATEMENT  CQNG
  IF ( Kprint>=2 ) WRITE (Lun,'(''1QNG QUICK CHECK''/)')
  !
  ! TEST ON IER = 0
  !
  Ipass = 1
  epsabs = 0.0E+00
  epmach = R1MACH(4)
  uflow = R1MACH(1)
  epsrel = MAX(SQRT(epmach),0.1E-07)
  a = 0.0E+00
  b = 0.1E+01
  CALL QNG(F1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
  ierv(1) = ier
  ip = 0
  error = ABS(exact1-result)
  IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
  IF ( ip==0 ) Ipass = 0
  IF ( Kprint/=0 ) CALL CPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,&
    ierv,1)
  !
  ! TEST ON IER = 1
  !
  CALL QNG(F2N,a,b,uflow,0.0E+00,result,abserr,neval,ier)
  ierv(1) = ier
  ip = 0
  IF ( ier==1 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  IF ( Kprint/=0 ) CALL CPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,&
    ierv,1)
  !
  ! TEST ON IER = 6
  !
  epsabs = 0.0E+00
  epsrel = 0.0E+00
  CALL QNG(F1N,a,b,epsabs,0.0E+00,result,abserr,neval,ier)
  ierv(1) = ier
  ip = 0
  IF ( ier==6.AND.result==0.0E+00.AND.abserr==0.0E+00.AND.neval==0 ) ip = 1
  IF ( ip==0 ) Ipass = 0
  IF ( Kprint/=0 ) CALL CPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,&
    ierv,1)
  !
  IF ( Kprint>=1 ) THEN
    IF ( Ipass==0 ) THEN
      WRITE (Lun,'(/'' SOME TEST(S) IN CQNG FAILED''/)')
    ELSEIF ( Kprint>=2 ) THEN
      WRITE (Lun,'(/'' ALL TEST(S) IN CQNG PASSED''/)')
    ENDIF
  ENDIF
END SUBROUTINE CQNG
