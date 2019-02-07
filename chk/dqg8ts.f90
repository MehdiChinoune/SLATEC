!*==DQG8TS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DQG8TS
SUBROUTINE DQG8TS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !*--DQG8TS5
  !***BEGIN PROLOGUE  DQG8TS
  !***PURPOSE  Quick check for DGAUS8.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QG8TST-S, DQG8TS-D)
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  D1MACH, DFQD1, DFQD2, DGAUS8, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
  !   910501  Added PURPOSE and TYPE records.  (WRB)
  !   910708  Minor modifications in use of KPRINT.  (WRB)
  !   920213  Code restructured to test DGAUS8 for all values of KPRINT,
  !           second accuracy test added and testing of error returns
  !           revised.  (WRB)
  !***END PROLOGUE  DQG8TS
  !     .. Scalar Arguments ..
  INTEGER Ipass , Kprint , Lun
  !     .. Local Scalars ..
  INTEGER ierr , kontrl
  REAL(8) :: a , ans , b , cor , err , req , tol
  LOGICAL fatal
  !     .. External Functions ..
  REAL(8) :: D1MACH , DFQD1 , DFQD2
  EXTERNAL D1MACH , DFQD1 , DFQD2
  !     .. External Subroutines ..
  EXTERNAL DGAUS8 , XGETF , XSETF
  !     .. Intrinsic Functions ..
  INTRINSIC ABS , ATAN , EXP , SQRT
  !***FIRST EXECUTABLE STATEMENT  DQG8TS
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
  CALL DGAUS8(DFQD1,a,b,err,ans,ierr)
  cor = 2.0D0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED' , a , b , ans , cor , err , &
      ierr
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED' , a , b , ans , cor , err , &
      ierr
  ENDIF
  !
  !     Second accuracy test.
  !
  a = 0.0D0
  b = 4.0D0*ATAN(1.0D0)
  err = tol/100.0D0
  CALL DGAUS8(DFQD2,a,b,err,ans,ierr)
  cor = (EXP(b)-1.0D0)/101.0D0
  IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
    IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED' , a , b , ans , cor , err , &
      ierr
  ELSE
    Ipass = 0
    IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED' , a , b , ans , cor , err , &
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
  a = 0.0D0
  b = 1.0D0
  cor = 2.0D0
  err = 100.0D0*D1MACH(4)
  req = err
  CALL DGAUS8(DFQD1,a,b,err,ans,ierr)
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
  !     Test DGAUS8 with A and B nearly equal.
  !
  a = 2.0D0
  b = a*(1.0D0+D1MACH(4))
  cor = 0.0D0
  err = tol
  !
  CALL DGAUS8(DFQD1,a,b,err,ans,ierr)
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
  99003 FORMAT ('1'/' DGAUS8 Quick Check')
  99004 FORMAT (/' Accuracy test of DGAUS8 ',A/' A = ',F10.5,'   B = ',&
    F10.5/' Computed result = ',D14.7,'   Exact result = ',&
    D14.7/' Tolerance = ',D14.7,'   IERR = ',I2/)
  99005 FORMAT (/' Test error returns'/' 2 error messages expected'/)
  99006 FORMAT (' Test of DGAUS8 ',A/' REQ =',D10.2,5X,'ANS =',D20.13,5X,'IERR =',&
    I2,5X,'should be 2'/' ERR =',D10.2,' CORRECT =',D20.13/)
  99007 FORMAT (' Test of A and B nearly equal ',A)
  99008 FORMAT (/,' ***************DGAUS8 PASSED ALL TESTS**************')
  99009 FORMAT (/,' ***************DGAUS8 FAILED SOME TESTS*************')
END SUBROUTINE DQG8TS
