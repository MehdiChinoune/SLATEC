MODULE TEST42_MOD
  IMPLICIT NONE

CONTAINS
  !** DAVNTS
  SUBROUTINE DAVNTS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DAVINT.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (AVNTST-S, DAVNTS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DAVINT, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920210  Code restructured and revised to test error returns for all
    !           values of KPRINT.  (WRB)

    INTEGER kontrl
    REAL(8) :: D1MACH
    INTEGER i, ierr, Ipass, Kprint, Lun, n
    REAL(8) :: a, ans, b, del, rn1, sqb, tol, tol1, x(501), xint, y(501)
    LOGICAL fatal
    !* FIRST EXECUTABLE STATEMENT  DAVNTS
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/' DAVINT Quick Check')
    Ipass = 1
    tol = MAX(.0001D0,SQRT(D1MACH(4)))
    tol1 = 1.0D-2*tol
    !
    !     Perform first accuracy test.
    !
    a = 0.0D0
    b = 5.0D0
    xint = EXP(5.0D0) - 1.0D0
    n = 500
    rn1 = n - 1
    sqb = SQRT(b)
    del = 0.4D0*(b-a)/(n-1)
    DO i = 1, n
      x(i) = sqb*SQRT(a+(i-1)*(b-a)/rn1) + del
      y(i) = EXP(x(i))
    END DO
    CALL DAVINT(x,y,n,a,b,ans,ierr)
    !
    !     See if test was passed.
    !
    IF ( ABS(ans-xint)>tol ) THEN
      Ipass = 0
      IF ( Kprint>=3 ) WRITE (Lun,99009) ierr, ans, xint
    END IF
    !
    !     Perform second accuracy test.
    !
    x(1) = 0.0D0
    x(2) = 5.0D0
    y(1) = 1.0D0
    y(2) = 0.5D0
    a = -0.5D0
    b = 0.5D0
    xint = 1.0D0
    CALL DAVINT(x,y,2,a,b,ans,ierr)
    !
    !     See if test was passed.
    !
    IF ( ABS(ans-xint)>tol1 ) THEN
      Ipass = 0
      IF ( Kprint>=3 ) WRITE (Lun,99009) ierr, ans, xint
    END IF
    !
    !     Send message indicating passage or failure of tests.
    !
    IF ( Kprint>=2 ) THEN
      IF ( Ipass==1 ) THEN
        IF ( Kprint>=3 ) WRITE (Lun,99002)
        99002 FORMAT (/' DAVINT passed both accuracy tests.')
      ELSE
        WRITE (Lun,99003)
        99003 FORMAT (/' DAVINT failed at least one accuracy test.')
      END IF
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    fatal = .FALSE.
    CALL XERCLR
    !
    IF ( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (/' Test error returns from DAVINT'/' 4 error messages expected'/ )
    END IF
    DO i = 1, 20
      x(i) = (i-1)/19.0D0 - 0.01D0
      IF ( i/=1 ) y(i) = x(i)/(EXP(x(i))-1.0)
    END DO
    !
    !     Test IERR = 1 error return.
    !
    y(1) = 1.0D0
    CALL DAVINT(x,y,20,0.0D0,1.0D0,ans,ierr)
    IF ( ierr/=1 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99010) ierr, 1
    END IF
    CALL XERCLR
    !
    !     Test IERR = 2 error return.
    !
    CALL DAVINT(x,y,20,1.0D0,0.0D0,ans,ierr)
    IF ( ierr/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99010) ierr, 2
    END IF
    IF ( ans/=0.0D0 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 5 error return.
    !
    CALL DAVINT(x,y,1,0.0D0,1.0D0,ans,ierr)
    IF ( ierr/=5 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99010) ierr, 5
    END IF
    IF ( ans/=0.0D0 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 4 error return.
    !
    x(1) = 1.0D0/19.0D0
    x(2) = 0.0D0
    CALL DAVINT(x,y,20,0.0D0,1.0D0,ans,ierr)
    IF ( ierr/=4 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99010) ierr, 4
    END IF
    IF ( ans/=0.0D0 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 3 error return.
    !
    x(1) = 0.0D0
    x(2) = 1.0D0/19.0D0
    CALL DAVINT(x,y,20,0.0D0,.01D0,ans,ierr)
    IF ( ierr/=3 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99010) ierr, 3
    END IF
    IF ( ans/=0.0D0 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF ( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Reset XERMSG control variables and write summary.
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99005)
        99005 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99006)
      99006 FORMAT (/' All incorrect argument tests PASSED')
    END IF
    !
    !     Write PASS/FAIL message.
    !
    IF ( Ipass==1.AND.Kprint>=3 ) WRITE (Lun,99007)
    99007 FORMAT (/' ***************DAVINT PASSED ALL TESTS***************')
    IF ( Ipass==0.AND.Kprint>=2 ) WRITE (Lun,99008)
    99008 FORMAT (/' ***************DAVINT FAILED SOME TESTS**************')
    RETURN
    99009 FORMAT (/' FAILED ACCURACY TEST'/' IERR=',I2,5X,'COMPUTED ANS=',&
      E20.11/14X,'CORRECT ANS=',D20.11,5X,'REQUESTED ERR=',D10.2)
    99010 FORMAT (/' IERR =',I2,' and it should =',I2/)
    99011 FORMAT (1X,'ANS .NE. 0')
  END SUBROUTINE DAVNTS
  !** DQG8TS
  SUBROUTINE DQG8TS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DGAUS8.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QG8TST-S, DQG8TS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DFQD1, DFQD2, DGAUS8, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920213  Code restructured to test DGAUS8 for all values of KPRINT,
    !           second accuracy test added and testing of error returns
    !           revised.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER ierr, kontrl
    REAL(8) :: a, ans, b, cor, err, req, tol
    LOGICAL fatal
    !     .. External Functions ..
    REAL(8), EXTERNAL :: D1MACH
    !     .. External Subroutines ..
    EXTERNAL :: DGAUS8, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, ATAN, EXP, SQRT
    !* FIRST EXECUTABLE STATEMENT  DQG8TS
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
      IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr
    END IF
    !
    !     Second accuracy test.
    !
    a = 0.0D0
    b = 4.0D0*ATAN(1.0D0)
    err = tol/100.0D0
    CALL DGAUS8(DFQD2,a,b,err,ans,ierr)
    cor = (EXP(b)-1.0D0)/101.0D0
    IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
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
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, err, cor
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99006) 'FAILED', req, ans, ierr, err, cor
      Ipass = 0
      fatal = .TRUE.
    END IF
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
    END IF
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99001)
        99001 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/' All incorrect argument tests PASSED')
    END IF
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
  !** DQN79Q
  SUBROUTINE DQN79Q(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for DQNC79.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QN79QX-S, DQN79Q-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  D1MACH, DFQD1, DFQD2, DQNC79, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920213  Code restructured to test DQNC79 for all values of KPRINT,
    !           second accuracy test added and testing of error returns
    !           revised.  (WRB)

    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER ierr, kontrl, nfct
    REAL(8) :: a, ans, b, cor, err, req, tol
    LOGICAL fatal
    !     .. External Functions ..
    REAL(8), EXTERNAL :: D1MACH
    !     .. External Subroutines ..
    EXTERNAL :: DQNC79, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, SQRT
    !* FIRST EXECUTABLE STATEMENT  DQN79Q
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
      IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr, nfct
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr, nfct
    END IF
    !
    !     Second accuracy test.
    !
    a = 0.0D0
    b = 4.0D0*ATAN(1.0D0)
    err = tol/10.0D0
    CALL DQNC79(DFQD2,a,b,err,ans,ierr,nfct)
    cor = (EXP(b)-1.0D0)/101.0D0
    IF ( ABS(ans-cor)<=tol.AND.ierr==1 ) THEN
      IF ( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr, nfct
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr, nfct
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF ( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
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
      IF ( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, err, cor
    ELSE
      IF ( Kprint>=2 ) WRITE (Lun,FMT=99006) 'FAILED', req, ans, ierr, err, cor
      Ipass = 0
      fatal = .TRUE.
    END IF
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
    END IF
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99001)
        99001 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/' All incorrect argument tests PASSED')
    END IF
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
  !** DFQD1
  REAL(8) FUNCTION DFQD1(X)
    IMPLICIT NONE
    !>
    !***
    !  Function evaluator for DQNC79 and DGAUS8 quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FQD1-S, DFQD1-D)
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !***
    ! **See also:**  DQG8TS, DQN79Q
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   920229  DATE WRITTEN

    !     .. Scalar Arguments ..
    REAL(8), INTENT(IN) :: X
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  DFQD1
    DFQD1 = 0.0D0
    IF ( X>0.0D0 ) DFQD1 = 1.0D0/SQRT(X)
  END FUNCTION DFQD1
  !** DFQD2
  REAL(8) FUNCTION DFQD2(X)
    IMPLICIT NONE
    !>
    !***
    !  Function evaluator for DQNC79 and DGAUS8 quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (FQD2-S, DFQD2-D)
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !***
    ! **See also:**  DQG8TS, DQN79Q
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   920229  DATE WRITTEN

    !     .. Scalar Arguments ..
    REAL(8), INTENT(IN) :: X
    !     .. Intrinsic Functions ..
    INTRINSIC COS, EXP
    !* FIRST EXECUTABLE STATEMENT  DFQD2
    DFQD2 = EXP(X)*COS(10.0D0*X)
  END FUNCTION DFQD2
END MODULE TEST42_MOD
!** TEST42
PROGRAM TEST42
  USE TEST42_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2
  !***
  ! **Type:**      DOUBLE PRECISION (TEST41-S, TEST42-D)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        DAVINT   DGAUS8   DQNC79
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DAVNTS, DQG8TS, DQN79Q, I1MACH, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST42
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test DAVINT
  !
  CALL DAVNTS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DGAUS8
  !
  CALL DQG8TS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test DQNC79
  !
  CALL DQN79Q(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST42 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST42 *************')
  END IF
  STOP
END PROGRAM TEST42
