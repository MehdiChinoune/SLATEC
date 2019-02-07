!*==FZTEST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FZTEST
SUBROUTINE FZTEST(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--FZTEST5
  !***BEGIN PROLOGUE  FZTEST
  !***PURPOSE  Quick check for FZERO.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FZTEST-S, DFZTST-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  FZERO, R1MACH, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920212  Code completely restructured to test IFLAG for all values
  !           of KPRINT.  (WRB)
  !***END PROLOGUE  FZTEST
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Local Scalars ..
  INTEGER iflag , kontrl
  REAL ae , b , c , pi , r , re , tol
  LOGICAL fatal
  !     .. External Functions ..
  REAL R1MACH
  EXTERNAL R1MACH
  !     .. External Subroutines ..
  EXTERNAL FZERO , XERCLR , XGETF , XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , ATAN , MAX , SIN , SQRT , TAN
  !***FIRST EXECUTABLE STATEMENT  FZTEST
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1'/' FZERO QUICK CHECK')
  Ipass = 1
  pi = 4.0E0*ATAN(1.0E0)
  re = 1.0E-6
  ae = 1.0E-6
  tol = MAX(1.0E-5,SQRT(R1MACH(4)))
  !
  !     Set up and solve example problem
  !
  b = 0.1E0
  c = 4.0E0
  r = c - b
  CALL FZERO(SIN,b,c,r,re,ae,iflag)
  !
  !     See if test was passed.
  !
  IF ( ABS(b-pi)<=tol.AND.ABS(c-pi)<=tol ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED' , b , c , iflag
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED' , b , c , iflag
  ENDIF
  !
  !     Trigger 2 error conditions
  !
  CALL XGETF(kontrl)
  IF ( Kprint<=2 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  fatal = .FALSE.
  CALL XERCLR
  !
  IF ( Kprint>=3 ) WRITE (Lun,99002)
  99002 FORMAT (/' IFLAG 3 and 4 tests')
  b = 1.0E0
  !
  !     IFLAG=3 (Singular point)
  !
  c = 2.0E0
  r = 0.5E0*(b+c)
  CALL FZERO(TAN,b,c,b,re,ae,iflag)
  IF ( iflag/=3 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99008) iflag , 2
  ENDIF
  !
  !     IFLAG=4 (No sign change)
  !
  b = -3.0E0
  c = -0.1E0
  r = 0.5E0*(b+c)
  CALL FZERO(SIN,b,c,r,re,ae,iflag)
  IF ( iflag/=4 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99008) iflag , 4
  ENDIF
  !
  CALL XERCLR
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99003)
      99003     FORMAT (/' At least IFLAG test failed')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99004)
    99004   FORMAT (/' All IFLAG tests passed')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99005)
  99005 FORMAT (/' ***************FZERO PASSED ALL TESTS**************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
  99006 FORMAT (/' ***************FZERO FAILED SOME TESTS*************')
  RETURN
  99007 FORMAT (' Accuracy test ',&
    A/' Example problem results:  (answer = PI),  B =',F20.14,' C =',&
    F20.14/' IFLAG =',I2)
  99008 FORMAT (/' IFLAG test FAILED.  IFLAG =',I2,', but should ','have been',I2)
END SUBROUTINE FZTEST
