!*==DFZTST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFZTST
SUBROUTINE DFZTST(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DFZTST5
  !***BEGIN PROLOGUE  DFZTST
  !***PURPOSE  Quick check for DFZERO.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FZTEST-S, DFZTST-D)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***ROUTINES CALLED  D1MACH, DFZERO, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   920212  DATE WRITTEN
  !***END PROLOGUE  DFZTST
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  INTEGER iflag, kontrl
  REAL(8) :: ae, b, c, pi, r, re, tol
  LOGICAL fatal
  !     .. External Functions ..
  REAL(8) :: D1MACH
  EXTERNAL D1MACH
  !     .. External Subroutines ..
  EXTERNAL DFZERO, XERCLR, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, ATAN, DSIN, DTAN, MAX, SQRT
  !***FIRST EXECUTABLE STATEMENT  DFZTST
  IF ( Kprint>=2 ) WRITE (Lun,99001)
  99001 FORMAT ('1'/' DFZERO QUICK CHECK')
  Ipass = 1
  pi = 4.0D0*ATAN(1.0D0)
  re = 1.0D-10
  ae = 1.0D-10
  tol = MAX(1.0D-9,SQRT(D1MACH(4)))
  !
  !     Set up and solve example problem
  !
  b = 0.1D0
  c = 4.0D0
  r = c - b
  CALL DFZERO(DSIN,b,c,r,re,ae,iflag)
  !
  !     See if test was passed.
  !
  IF ( ABS(b-pi)<=tol.AND.ABS(c-pi)<=tol ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED', b, c, iflag
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED', b, c, iflag
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
  b = 1.0D0
  !
  !     IFLAG=3 (Singular point)
  !
  c = 2.0D0
  r = 0.5D0*(b+c)
  CALL DFZERO(DTAN,b,c,b,re,ae,iflag)
  IF ( iflag/=3 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 2
  ENDIF
  !
  !     IFLAG=4 (No sign change)
  !
  b = -3.0D0
  c = -0.1D0
  r = 0.5D0*(b+c)
  CALL DFZERO(DSIN,b,c,r,re,ae,iflag)
  IF ( iflag/=4 ) THEN
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99008) iflag, 4
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
  99005 FORMAT (/' ***************DFZERO PASSED ALL TESTS**************')
  IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99006)
  99006 FORMAT (/' ***************DFZERO FAILED SOME TESTS*************')
  RETURN
  99007 FORMAT (' Accuracy test ',&
    A/' Example problem results:  (answer = PI),  B =',F20.14,' C =',&
    F20.14/' IFLAG =',I2)
  99008 FORMAT (/' IFLAG test FAILED.  IFLAG =',I2,', but should ','have been',I2)
END SUBROUTINE DFZTST
