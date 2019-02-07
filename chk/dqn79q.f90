!*==DQN79Q.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DQN79Q
SUBROUTINE DQN79Q(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DQN79Q5
  !***BEGIN PROLOGUE  DQN79Q
  !***PURPOSE  Quick check for DQNC79.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QN79QX-S, DQN79Q-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DFQD1, DFQD2, DQNC79, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920213  Code restructured to test DQNC79 for all values of KPRINT,
  !           second accuracy test added and testing of error returns
  !           revised.  (WRB)
  !***END PROLOGUE  DQN79Q
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Local Scalars ..
  INTEGER ierr , kontrl , nfct
  DOUBLE PRECISION a , ans , b , cor , err , req , tol
  LOGICAL fatal
  !     .. External Functions ..
  DOUBLE PRECISION D1MACH , DFQD1 , DFQD2
  EXTERNAL D1MACH , DFQD1 , DFQD2
  !     .. External Subroutines ..
  EXTERNAL DQNC79 , XGETF , XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , MAX , SQRT
  !***FIRST EXECUTABLE STATEMENT  DQN79Q
  IF ( Kprint>=2 ) WRITE (Lun,FMT=99003)
  !
  !     Initialize variables for testing.
  !
  tol = SQRT(D1MACH(4))
  Ipass = 1
  !
  !     First accuracy test.
  !
  a = 1.0D0
  b = 4.0D0
  err = tol/100.0D0
  CALL DQNC79(DFQD1,a,b,err,ans,ierr,nfct)
  cor = 2.0D0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED' , a , b , ans , cor , err , &
      ierr , nfct
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED' , a , b , ans , cor , err , &
      ierr , nfct
  ENDIF
  !
  !     Second accuracy test.
  !
  a = 0.0D0
  b = 4.0D0*ATAN(1.0D0)
  err = tol/10.0D0
  CALL DQNC79(DFQD2,a,b,err,ans,ierr,nfct)
  cor = (EXP(b)-1.0D0)/101.0D0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED' , a , b , ans , cor , err , &
      ierr , nfct
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED' , a , b , ans , cor , err , &
      ierr , nfct
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
  a = 0.0D0
  b = 1.0D0
  cor = 2.0D0
  err = 100.0D0*D1MACH(4)
  req = err
  CALL DQNC79(DFQD1,a,b,err,ans,ierr,nfct)
  !
  !     See if test passed.
  !
  IF ( ierr==2 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED' , req , ans , ierr , &
      err , cor
  ELSE
    IF ( Kprint>=2 ) WRITE (Lun,FMT=99006) 'FAILED' , req , ans , ierr , &
      err , cor
    Ipass = 0
    fatal = .TRUE.
  ENDIF
  !
  !     Test DQNC79 with A and B nearly equal.
  !
  a = 2.0D0
  b = a*(1.0D0+D1MACH(4))
  cor = 0.0D0
  err = tol
  !
  CALL DQNC79(DFQD1,a,b,err,ans,ierr,nfct)
  !
  !     Check to see if test passed.
  !
  IF ( ierr==-1.AND.ans==0.0D0 ) THEN
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
  99003 FORMAT ('1'/' DQNC79 Quick Check')
  99004 FORMAT (/' Accuracy test of DQNC79 ',A/' A = ',F10.5,'   B = ',&
    F10.5/' Computed result = ',D14.7,'   Exact result = ',&
    D14.7/' Tolerance = ',D14.7,'   IERR = ',I2,&
    '   Number of function evals = ',I5/)
  99005 FORMAT (/' Test error returns'/' 2 error messages expected'/)
  99006 FORMAT (' Test of DQNC79 ',A/' REQ =',D10.2,5X,'ANS =',D20.13,5X,'IERR =',&
    I2,5X,'should be 2'/' ERR =',D10.2,' CORRECT =',D20.13/)
  99007 FORMAT (' Test of A and B nearly equal ',A)
  99008 FORMAT (/' ***************DQNC79 PASSED ALL TESTS***************')
  99009 FORMAT (/' ***************DQNC79 FAILED SOME TESTS**************')
END SUBROUTINE DQN79Q
