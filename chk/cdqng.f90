!*==CDQNG.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CDQNG
      SUBROUTINE CDQNG(Lun,Kprint,Ipass)
      IMPLICIT NONE
!*--CDQNG5
!*** Start of declarations inserted by SPAG
      INTEGER Lun
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CDQNG
!***PURPOSE  Quick check for DQNG.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (CQNG-S, CDQNG-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  D1MACH, DF1N, DF2N, DPRIN, DQNG
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901205  Added PASS/FAIL message and changed the name of the first
!           argument.  (RWC)
!   910501  Added PURPOSE and TYPE records.  (WRB)
!***END PROLOGUE  CDQNG
!
! FOR FURTHER DOCUMENTATION SEE ROUTINE CQPDOC
!
      DOUBLE PRECISION a , abserr , b , D1MACH , epmach , epsabs , epsrel , 
     &                 exact1 , error , exact2 , DF1N , DF2N , result , uflow
      INTEGER ier , ierv , ip , Ipass , Kprint , neval
      DIMENSION ierv(1)
      EXTERNAL DF1N , DF2N
      DATA exact1/0.7281029132255818D+00/
      DATA exact2/0.1D+02/
!***FIRST EXECUTABLE STATEMENT  CDQNG
      IF ( Kprint>=2 ) WRITE (Lun,'(''1DQNG QUICK CHECK''/)')
!
! TEST ON IER = 0
!
      Ipass = 1
      epsabs = 0.0D+00
      epmach = D1MACH(4)
      uflow = D1MACH(1)
      epsrel = MAX(SQRT(epmach),0.1D-07)
      a = 0.0D+00
      b = 0.1D+01
      CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
      CALL DQNG(DF1N,a,b,epsabs,epsrel,result,abserr,neval,ier)
      ierv(1) = ier
      ip = 0
      error = ABS(exact1-result)
      IF ( ier==0.AND.error<=abserr.AND.abserr<=epsrel*ABS(exact1) ) ip = 1
      IF ( ip==0 ) Ipass = 0
      IF ( Kprint/=0 ) CALL DPRIN(Lun,0,Kprint,ip,exact1,result,abserr,neval,
     &                            ierv,1)
!
! TEST ON IER = 1
!
      CALL DQNG(DF2N,a,b,uflow,0.0D+00,result,abserr,neval,ier)
      ierv(1) = ier
      ip = 0
      IF ( ier==1 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      IF ( Kprint/=0 ) CALL DPRIN(Lun,1,Kprint,ip,exact2,result,abserr,neval,
     &                            ierv,1)
!
! TEST ON IER = 6
!
      epsabs = 0.0D+00
      epsrel = 0.0D+00
      CALL DQNG(DF1N,a,b,epsabs,0.0D+00,result,abserr,neval,ier)
      ierv(1) = ier
      ip = 0
      IF ( ier==6.AND.result==0.0D+00.AND.abserr==0.0D+00.AND.neval==0 ) ip = 1
      IF ( ip==0 ) Ipass = 0
      IF ( Kprint/=0 ) CALL DPRIN(Lun,6,Kprint,ip,exact1,result,abserr,neval,
     &                            ierv,1)
!
      IF ( Kprint>=1 ) THEN
        IF ( Ipass==0 ) THEN
          WRITE (Lun,'(/'' SOME TEST(S) IN CDQNG FAILED''/)')
        ELSEIF ( Kprint>=2 ) THEN
          WRITE (Lun,'(/'' ALL TEST(S) IN CDQNG PASSED''/)')
        ENDIF
      ENDIF
      END SUBROUTINE CDQNG
