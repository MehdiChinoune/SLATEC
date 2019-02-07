!*==QG8TST.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK QG8TST
SUBROUTINE QG8TST(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--QG8TST5
  !***BEGIN PROLOGUE  QG8TST
  !***PURPOSE  Quick check for GAUS8.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QG8TST-S, DQG8TS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  FQD1, FQD2, GAUS8, R1MACH, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920213  Code restructured to test GAUS8 for all values of KPRINT,
  !           second accuracy test added and testing of error returns
  !           revised.  (WRB)
  !***END PROLOGUE  QG8TST
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  INTEGER ierr, kontrl
  REAL a, ans, b, cor, err, req, tol
  LOGICAL fatal
  !     .. External Functions ..
  REAL FQD1, FQD2, R1MACH
  EXTERNAL FQD1, FQD2, R1MACH
  !     .. External Subroutines ..
  EXTERNAL GAUS8, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, ATAN, EXP, SQRT
  !***FIRST EXECUTABLE STATEMENT  QG8TST
  IF ( Kprint>=2 ) WRITE (Lun,FMT=99003)
  !
  !     Initialize variables for testing.
  !
  tol = SQRT(R1MACH(4))
  Ipass = 1
  !
  !     First accuracy test.
  !
  a = 1.0E0
  b = 4.0E0
  err = tol/100.0E0
  CALL GAUS8(FQD1,a,b,err,ans,ierr)
  cor = 2.0E0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, &
      ierr
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, &
      ierr
  ENDIF
  !
  !     Second accuracy test.
  !
  a = 0.0E0
  b = 4.0E0*ATAN(1.0E0)
  err = tol/100.0E0
  CALL GAUS8(FQD2,a,b,err,ans,ierr)
  cor = (EXP(b)-1.0E0)/101.0E0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, &
      ierr
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, &
      ierr
  ENDIF
  !
  !     Test error returns.
  !
  CALL XGETF(kontrl)
  IF ( Kprint<=2 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  fatal = .FALSE.
  !
  IF ( Kprint>=3 ) WRITE (Lun,FMT=99005)
  !
  !     Test with a discontinuous integrand and a tight error tolerance.
  !
  a = 0.0E0
  b = 1.0E0
  cor = 2.0E0
  err = 100.0E0*R1MACH(4)
  req = err
  CALL GAUS8(FQD1,a,b,err,ans,ierr)
  !
  !     See if test passed.
  !
  IF ( ierr==2 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, &
      err, cor
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, &
      err, cor
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  !
  !     Test GAUS8 with A and B nearly equal.
  !
  a = 2.0E0
  b = a*(1.0E0+R1MACH(4))
  cor = 0.0E0
  err = tol
  !
  CALL GAUS8(FQD1,a,b,err,ans,ierr)
  !
  !     Check to see if test passed.
  !
  IF ( ierr==-1.AND.ans==0.0E0 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99007) 'PASSED'
  ELSE
    Ipass = 0
    fatal = .TRUE.
    IF ( Kprint>=2 ) WRITE (Lun,99007) 'FAILED'
  ENDIF
  !
  CALL XSETF(kontrl)
  IF ( fatal ) THEN
    IF ( Kprint>=2 ) THEN
      WRITE (Lun,99001)
      99001     FORMAT (/' At least one incorrect argument test FAILED')
    ENDIF
  ELSEIF ( Kprint>=3 ) THEN
    WRITE (Lun,99002)
    99002   FORMAT (/' All incorrect argument tests PASSED')
  ENDIF
  !
  IF ( Ipass==1.AND.Kprint>=3 ) WRITE (Lun,FMT=99008)
  IF ( Ipass==0.AND.Kprint>=2 ) WRITE (Lun,FMT=99009)
  RETURN
  !
  99003 FORMAT ('1'/' GAUS8 Quick Check')
  99004 FORMAT (/' Accuracy test of GAUS8 ',A/' A = ',F10.5,'   B = ',&
    F10.5/' Computed result = ',E14.7,'   Exact result = ',&
    E14.7/' Tolerance = ',E14.7,'   IERR = ',I2/)
  99005 FORMAT (/' Test error returns'/' 2 error messages expected'/)
  99006 FORMAT (' Test of GAUS8 ',A/' REQ =',E10.2,5X,'ANS =',E20.13,5X,'IERR =',&
    I2,5X,'should be 2'/' ERR =',E10.2,' CORRECT =',E20.13/)
  99007 FORMAT (' Test of A and B nearly equal ',A)
  99008 FORMAT (/,' ***************GAUS8 PASSED ALL TESTS***************')
  99009 FORMAT (/,' ***************GAUS8 FAILED SOME TESTS**************')
END SUBROUTINE QG8TST
