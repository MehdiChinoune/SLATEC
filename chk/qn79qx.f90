!DECK QN79QX
SUBROUTINE QN79QX(Lun,Kprint,Ipass)
  IMPLICIT NONE
  INTEGER kontrl
  !***BEGIN PROLOGUE  QN79QX
  !***PURPOSE  Quick check for QNC79.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QN79QX-S, DQN79Q-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  FQD1, FQD2, QNC79, R1MACH, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920213  Code restructured to test QNC79 for all values of KPRINT,
  !           second accuracy test added and testing of error returns
  !           revised.  (WRB)
  !***END PROLOGUE  QN79QX
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  INTEGER ierr, nfct
  REAL a, ans, b, cor, err, req, tol
  LOGICAL fatal
  !     .. External Functions ..
  REAL FQD1, FQD2, R1MACH
  EXTERNAL FQD1, FQD2, R1MACH
  !     .. External Subroutines ..
  EXTERNAL QNC79, XGETF, XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, MAX, SQRT
  !***FIRST EXECUTABLE STATEMENT  QN79QX
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
  CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
  cor = 2.0E0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, &
      ierr, nfct
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, &
      ierr, nfct
  ENDIF
  !
  !     Second accuracy test.
  !
  a = 0.0E0
  b = 4.0E0*ATAN(1.0E0)
  err = tol/10.0E0
  CALL QNC79(FQD2,a,b,err,ans,ierr,nfct)
  cor = (EXP(b)-1.0E0)/101.0E0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, &
      ierr, nfct
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, &
      ierr, nfct
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
  CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
  !
  !     See if test passed.
  !
  IF ( ierr==2 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, &
      err, cor
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99006) 'FAILED', req, ans, ierr, &
      err, cor
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  !
  !     Test QNC79 with A and B nearly equal.
  !
  a = 2.0E0
  b = a*(1.0E0+R1MACH(4))
  cor = 0.0E0
  err = tol
  !
  CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
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
  99003 FORMAT ('1'/' QNC79 Quick Check')
  99004 FORMAT (/' Accuracy test of QNC79 ',A/' A = ',F10.5,'   B = ',&
    F10.5/' Computed result = ',E14.7,'   Exact result = ',&
    E14.7/' Tolerance = ',E14.7,'   IERR = ',I2,&
    '   Number of function evals = ',I5/)
  99005 FORMAT (/' Test error returns'/' 2 error messages expected'/)
  99006 FORMAT (' Test of QNC79 ',A/' REQ =',E10.2,5X,'ANS =',E20.13,5X,'IERR =',&
    I2,5X,'should be 2'/' ERR =',E10.2,' CORRECT =',E20.13/)
  99007 FORMAT (' Test of A and B nearly equal ',A)
  99008 FORMAT (/' ***************QNC79 PASSED ALL TESTS****************')
  99009 FORMAT (/' ***************QNC79 FAILED SOME TESTS***************')
END SUBROUTINE QN79QX
