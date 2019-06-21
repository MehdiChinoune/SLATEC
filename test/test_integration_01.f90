MODULE TEST41_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** AVNTST
  SUBROUTINE AVNTST(Lun,Kprint,Ipass)
    !> Quick check for AVINT.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (AVNTST-S, DAVNTS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  AVINT, R1MACH, XERCLR, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   890911  Removed unnecessary intrinsics.  (WRB)
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920210  Code restructured and revised to test error returns for all
    !           values of KPRINT.  (WRB)
    USE slatec, ONLY : AVINT, R1MACH, XERCLR, XGETF, XSETF
    REAL(SP) :: a, ans, b, del, rn1, sqb, tol, tol1, x(501), xint, y(501)
    INTEGER :: i, ierr, Ipass, kontrl, Kprint, Lun, n
    LOGICAL :: fatal
    !* FIRST EXECUTABLE STATEMENT  AVNTST
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT ('1'/' AVINT Quick Check')
    Ipass = 1
    tol = MAX(.0001_SP,SQRT(R1MACH(4)))
    tol1 = 1.E-2_SP*tol
    !
    !     Perform first accuracy test.
    !
    a = 0._SP
    b = 5._SP
    xint = EXP(5._SP) - 1._SP
    n = 500
    rn1 = n - 1
    sqb = SQRT(b)
    del = 0.4_SP*(b-a)/(n-1)
    DO i = 1, n
      x(i) = sqb*SQRT(a+(i-1)*(b-a)/rn1) + del
      y(i) = EXP(x(i))
    END DO
    CALL AVINT(x,y,n,a,b,ans,ierr)
    !
    !     See if test was passed.
    !
    IF( ABS(ans-xint)>tol ) THEN
      Ipass = 0
      IF( Kprint>=3 ) WRITE (Lun,99009) ierr, ans, xint
    END IF
    !
    !     Perform second accuracy test.
    !
    x(1) = 0._SP
    x(2) = 5._SP
    y(1) = 1._SP
    y(2) = 0.5_SP
    a = -0.5_SP
    b = 0.5_SP
    xint = 1._SP
    CALL AVINT(x,y,2,a,b,ans,ierr)
    !
    !     See if test was passed.
    !
    IF( ABS(ans-xint)>tol1 ) THEN
      Ipass = 0
      IF( Kprint>=3 ) WRITE (Lun,99009) ierr, ans, xint
    END IF
    !
    !     Send message indicating passage or failure of tests.
    !
    IF( Kprint>=2 ) THEN
      IF( Ipass==1 ) THEN
        IF( Kprint>=3 ) WRITE (Lun,99002)
        99002 FORMAT (/' AVINT passed both accuracy tests.')
      ELSE
        WRITE (Lun,99003)
        99003 FORMAT (/' AVINT failed at least one accuracy test.')
      END IF
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    fatal = .FALSE.
    CALL XERCLR
    !
    IF( Kprint>=3 ) THEN
      WRITE (Lun,99004)
      99004 FORMAT (/' Test error returns from AVINT'/' 4 error messages expected'/)
    END IF
    DO i = 1, 20
      x(i) = (i-1)/19._SP - 0.01_SP
      IF( i/=1 ) y(i) = x(i)/(EXP(x(i))-1._SP)
    END DO
    !
    !     Test IERR = 1 error return.
    !
    y(1) = 1._SP
    CALL AVINT(x,y,20,0._SP,1._SP,ans,ierr)
    IF( ierr/=1 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99010) ierr, 1
    END IF
    CALL XERCLR
    !
    !     Test IERR = 2 error return.
    !
    CALL AVINT(x,y,20,1._SP,0._SP,ans,ierr)
    IF( ierr/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99010) ierr, 2
    END IF
    IF( ans/=0._SP ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 5 error return.
    !
    CALL AVINT(x,y,1,0._SP,1._SP,ans,ierr)
    IF( ierr/=5 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99010) ierr, 5
    END IF
    IF( ans/=0._SP ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 4 error return.
    !
    x(1) = 1._SP/19._SP
    x(2) = 0._SP
    CALL AVINT(x,y,20,0._SP,1._SP,ans,ierr)
    IF( ierr/=4 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99010) ierr, 4
    END IF
    IF( ans/=0._SP ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Test IERR = 3 error return.
    !
    x(1) = 0._SP
    x(2) = 1._SP/19._SP
    CALL AVINT(x,y,20,0._SP,.01_SP,ans,ierr)
    IF( ierr/=3 ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99010) ierr, 3
    END IF
    IF( ans/=0._SP ) THEN
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=3 ) WRITE (Lun,99011)
    END IF
    CALL XERCLR
    !
    !     Reset XERMSG control variables and write summary.
    !
    CALL XSETF(kontrl)
    IF( fatal ) THEN
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99005)
        99005 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99006)
      99006 FORMAT (/' All incorrect argument tests PASSED')
    END IF
    !
    !     Write PASS/FAIL message.
    !
    IF( Ipass==1 .AND. Kprint>=3 ) WRITE (Lun,99007)
    99007 FORMAT (/' ***************AVINT PASSED ALL TESTS***************')
    IF( Ipass==0 .AND. Kprint>=2 ) WRITE (Lun,99008)
    99008 FORMAT (/' ***************AVINT FAILED SOME TESTS**************')
    RETURN
    99009 FORMAT (/' FAILED ACCURACY TEST'/' IERR=',I2,5X,'COMPUTED ANS=',&
      E20.11/14X,'CORRECT ANS=',E20.11,5X,'REQUESTED ERR=',E10.2)
    99010 FORMAT (/' IERR =',I2,' and it should =',I2/)
    99011 FORMAT (1X,'ANS /= 0')
  END SUBROUTINE AVNTST
  !** QG8TST
  SUBROUTINE QG8TST(Lun,Kprint,Ipass)
    !> Quick check for GAUS8.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (QG8TST-S, DQG8TS-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  FQD1, FQD2, GAUS8, R1MACH, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920213  Code restructured to test GAUS8 for all values of KPRINT,
    !           second accuracy test added and testing of error returns
    !           revised.  (WRB)
    USE slatec, ONLY : GAUS8, R1MACH, XGETF, XSETF
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER :: ierr, kontrl
    REAL(SP) :: a, ans, b, cor, err, req, tol
    LOGICAL :: fatal
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, ATAN, EXP, SQRT
    !* FIRST EXECUTABLE STATEMENT  QG8TST
    IF( Kprint>=2 ) WRITE (Lun,FMT=99003)
    !
    !     Initialize variables for testing.
    !
    tol = SQRT(R1MACH(4))
    Ipass = 1
    !
    !     First accuracy test.
    !
    a = 1._SP
    b = 4._SP
    err = tol/100._SP
    CALL GAUS8(FQD1,a,b,err,ans,ierr)
    cor = 2._SP
    IF( ABS(ans-cor)<=tol .AND. ierr==1 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr
    END IF
    !
    !     Second accuracy test.
    !
    a = 0._SP
    b = 4._SP*ATAN(1._SP)
    err = tol/100._SP
    CALL GAUS8(FQD2,a,b,err,ans,ierr)
    cor = (EXP(b)-1._SP)/101._SP
    IF( ABS(ans-cor)<=tol .AND. ierr==1 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    fatal = .FALSE.
    !
    IF( Kprint>=3 ) WRITE (Lun,FMT=99005)
    !
    !     Test with a discontinuous integrand and a tight error tolerance.
    !
    a = 0._SP
    b = 1._SP
    cor = 2._SP
    err = 100._SP*R1MACH(4)
    req = err
    CALL GAUS8(FQD1,a,b,err,ans,ierr)
    !
    !     See if test passed.
    !
    IF( ierr==2 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, err, cor
    ELSE
      IF( Kprint>=2 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, err, cor
      Ipass = 0
      fatal = .TRUE.
    END IF
    !
    !     Test GAUS8 with A and B nearly equal.
    !
    a = 2._SP
    b = a*(1._SP+R1MACH(4))
    cor = 0._SP
    err = tol
    !
    CALL GAUS8(FQD1,a,b,err,ans,ierr)
    !
    !     Check to see if test passed.
    !
    IF( ierr==-1 .AND. ans==0._SP ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99007) 'PASSED'
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) WRITE (Lun,99007) 'FAILED'
    END IF
    !
    CALL XSETF(kontrl)
    IF( fatal ) THEN
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99001)
        99001 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/' All incorrect argument tests PASSED')
    END IF
    !
    IF( Ipass==1 .AND. Kprint>=3 ) WRITE (Lun,FMT=99008)
    IF( Ipass==0 .AND. Kprint>=2 ) WRITE (Lun,FMT=99009)
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
  !** QN79QX
  SUBROUTINE QN79QX(Lun,Kprint,Ipass)
    !> Quick check for QNC79.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (QN79QX-S, DQN79Q-D)
    !***
    ! **Author:**  (UNKNOWN)
    !***
    ! **Routines called:**  FQD1, FQD2, QNC79, R1MACH, XGETF, XSETF

    !* REVISION HISTORY  (YYMMDD)
    !   ??????  DATE WRITTEN
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
    !   910501  Added PURPOSE and TYPE records.  (WRB)
    !   910708  Minor modifications in use of KPRINT.  (WRB)
    !   920213  Code restructured to test QNC79 for all values of KPRINT,
    !           second accuracy test added and testing of error returns
    !           revised.  (WRB)
    USE slatec, ONLY : QNC79, R1MACH, XGETF, XSETF
    INTEGER :: kontrl
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    INTEGER :: ierr, nfct
    REAL(SP) :: a, ans, b, cor, err, req, tol
    LOGICAL :: fatal
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, MAX, SQRT
    !* FIRST EXECUTABLE STATEMENT  QN79QX
    IF( Kprint>=2 ) WRITE (Lun,FMT=99003)
    !
    !     Initialize variables for testing.
    !
    tol = SQRT(R1MACH(4))
    Ipass = 1
    !
    !     First accuracy test.
    !
    a = 1._SP
    b = 4._SP
    err = tol/100._SP
    CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
    cor = 2._SP
    IF( ABS(ans-cor)<=tol .AND. ierr==1 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr, nfct
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr, nfct
    END IF
    !
    !     Second accuracy test.
    !
    a = 0._SP
    b = 4._SP*ATAN(1._SP)
    err = tol/10._SP
    CALL QNC79(FQD2,a,b,err,ans,ierr,nfct)
    cor = (EXP(b)-1._SP)/101._SP
    IF( ABS(ans-cor)<=tol .AND. ierr==1 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99004) 'PASSED', a, b, ans, cor, err, ierr, nfct
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99004) 'FAILED', a, b, ans, cor, err, ierr, nfct
    END IF
    !
    !     Test error returns.
    !
    CALL XGETF(kontrl)
    IF( Kprint<=2 ) THEN
      CALL XSETF(0)
    ELSE
      CALL XSETF(1)
    END IF
    fatal = .FALSE.
    !
    IF( Kprint>=3 ) WRITE (Lun,FMT=99005)
    !
    !     Test with a discontinuous integrand and a tight error tolerance.
    !
    a = 0._SP
    b = 1._SP
    cor = 2._SP
    err = 100._SP*R1MACH(4)
    req = err
    CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
    !
    !     See if test passed.
    !
    IF( ierr==2 ) THEN
      IF( Kprint>=3 ) WRITE (Lun,FMT=99006) 'PASSED', req, ans, ierr, err, cor
    ELSE
      IF( Kprint>=2 ) WRITE (Lun,FMT=99006) 'FAILED', req, ans, ierr, err, cor
      Ipass = 0
      fatal = .TRUE.
    END IF
    !
    !     Test QNC79 with A and B nearly equal.
    !
    a = 2._SP
    b = a*(1._SP+R1MACH(4))
    cor = 0._SP
    err = tol
    !
    CALL QNC79(FQD1,a,b,err,ans,ierr,nfct)
    !
    !     Check to see if test passed.
    !
    IF( ierr==-1 .AND. ans==0._SP ) THEN
      IF( Kprint>=3 ) WRITE (Lun,99007) 'PASSED'
    ELSE
      Ipass = 0
      fatal = .TRUE.
      IF( Kprint>=2 ) WRITE (Lun,99007) 'FAILED'
    END IF
    !
    CALL XSETF(kontrl)
    IF( fatal ) THEN
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99001)
        99001 FORMAT (/' At least one incorrect argument test FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99002)
      99002 FORMAT (/' All incorrect argument tests PASSED')
    END IF
    !
    IF( Ipass==1 .AND. Kprint>=3 ) WRITE (Lun,FMT=99008)
    IF( Ipass==0 .AND. Kprint>=2 ) WRITE (Lun,FMT=99009)
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
  !** FQD1
  REAL(SP) FUNCTION FQD1(X)
    !> Function evaluator for QNC79 and GAUS8 quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (FQD1-S, DFQD1-D)
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !***
    ! **See also:**  QG8TST, QN79QX
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   920229  DATE WRITTEN

    !     .. Scalar Arguments ..
    REAL(SP),INTENT(IN) :: X
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !* FIRST EXECUTABLE STATEMENT  FQD1
    FQD1 = 0._SP
    IF( X>0._SP ) FQD1 = 1._SP/SQRT(X)
  END FUNCTION FQD1
  !** FQD2
  REAL(SP) FUNCTION FQD2(X)
    !> Function evaluator for QNC79 and GAUS8 quick checks.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (FQD2-S, DFQD2-D)
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !***
    ! **See also:**  QG8TST, QN79QX
    !***
    ! **Routines called:**  (NONE)

    !* REVISION HISTORY  (YYMMDD)
    !   920229  DATE WRITTEN

    !     .. Scalar Arguments ..
    REAL(SP),INTENT(IN) :: X
    !     .. Intrinsic Functions ..
    INTRINSIC COS, EXP
    !* FIRST EXECUTABLE STATEMENT  FQD2
    FQD2 = EXP(X)*COS(10._SP*X)
  END FUNCTION FQD2
END MODULE TEST41_MOD
!** TEST41
PROGRAM TEST41
  USE TEST41_MOD, ONLY : AVNTST, QG8TST, QN79QX
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  H2
  !***
  ! **Type:**      SINGLE PRECISION (TEST41-S, TEST42-D)
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
  !        AVINT    GAUS8    QNC79
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  AVNTST, I1MACH, QG8TST, QN79QX, XERMAX, XSETF,
  !                    XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST41
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test AVINT
  !
  CALL AVNTST(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test GAUS8
  !
  CALL QG8TST(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test QNC79
  !
  CALL QN79QX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST41 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST41 *************')
  END IF
  STOP
END PROGRAM TEST41
