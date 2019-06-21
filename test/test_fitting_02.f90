MODULE TEST28_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DLSEIT
  SUBROUTINE DLSEIT(Lun,Kprint,Ipass)
    !> Quick check for DLSEI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (LSEIQX-S, DLSEIT-D)
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Hanson, R. J., (SNLA)
    !           Haskell, Karen, (SNLA)
    !***
    ! **Description:**
    !
    !   The sample problem solved is from a paper by J. Stoer, in
    !   SIAM Journal of Numerical Analysis, June 1971.
    !
    !***
    ! **Routines called:**  D1MACH, DAXPY, DCOPY, DDOT, DLSEI, DNRM2, DVOUT

    !* REVISION HISTORY  (YYMMDD)
    !   790216  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
    !           to use D1MACH(4) rather than D1MACH(3) and cleaned up
    ! FORMATs.  (RWC)
    !   920722  Initialized IP(1) and IP(2) for CALL to DLSEI.  (BKS, WRB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    USE slatec, ONLY : D1MACH, DLSEI, DVOUT, XGETF, XSETF, XERCLR, NUMXER
    USE blas, ONLY : DAXPY
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(DP) :: cnorm, relerr(1), relnrm(1), resnrm(1), rnorme, rnorml(1), tnorm
    INTEGER :: i, idigit, jdigit, kontrl, ma, mdd, me, meap1, mep1, mg, &
      mode, n, nerr, np1
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(DP) :: d(11,6), err(5), prgopt(4), work(105), x(5)
    INTEGER :: ip(17)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !     .. Data statements ..
    !
    !     Define the data arrays for the example.  The array A contains
    !     the least squares equations.  (There are no equality constraints
    !     in this example).
    !
    REAL(DP), PARAMETER :: a(6,5) = RESHAPE( [ &
      -74._DP, 80._DP, 18._DP, -11._DP, -4._DP,   14._DP, -69._DP, 21._DP, 28._DP, 0._DP, &
      66._DP, -72._DP, -5._DP, 7._DP, 1._DP,     -12._DP, 66._DP, -30._DP, -23._DP, 3._DP, &
      3._DP, 8._DP, -7._DP, -4._DP, 1._DP,        4._DP, -12._DP, 4._DP, 4._DP, 0._DP ], &
      [6,5], ORDER = [2,1] )
    !
    !     The array G contains the inequality constraint equations,
    !     written in the sense
    !     (row vector)*(solution vector) >= (given value).
    !
    REAL(DP), PARAMETER :: g(5,5) = RESHAPE( [ -1._DP, -1._DP, -1._DP, -1._DP, -1._DP, &
      10._DP, 10._DP, -3._DP, 5._DP, 4._DP,    -8._DP, 1._DP, -2._DP, -5._DP, 3._DP, &
      8._DP, -1._DP, 2._DP, 5._DP, -3._DP,    -4._DP, -2._DP, 3._DP, -5._DP, 1._DP ], &
      [5,5], ORDER = [2,1] )
    !
    !     Define the least squares right-side vector.
    !
    REAL(DP), PARAMETER :: f(6) = [ -5._DP, -9._DP, 708._DP, 4165._DP, &
      -13266._DP, 8409._DP ]
    !
    !     Define the inequality constraint right-side vector.
    !
    REAL(DP), PARAMETER :: h(5) = [ -5._DP, 20._DP, -40._DP, 11._DP, -30._DP ]
    !
    !     Define the vector that is the known solution.
    !
    REAL(DP), PARAMETER :: sol(5) = [ 1._DP, 2._DP, -1._DP, 3._DP, -4._DP ]
    !* FIRST EXECUTABLE STATEMENT  DDLSEIT
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1TEST OF SUBROUTINE DLSEI')
    !
    !     Define the matrix dimensions, number of least squares equations,
    !     number of equality constraints, total number of equations, and
    !     number of variables.  Set ME=0 to indicate there are no equality
    !     constraints.
    !
    mdd = 11
    ma = 6
    mg = 5
    n = 5
    me = 0
    !
    ip(1) = 105
    ip(2) = 17
    !
    np1 = n + 1
    mep1 = me + 1
    meap1 = me + ma + 1
    !
    !     Copy the problem matrices.
    !
    d(ma+1:mdd,1:n) = g
    d(1:ma,1:n) = a
    !
    !     Copy the right-side vectors into the work array in compatible order.
    !
    d(ma+1:mdd,np1) = h
    d(1:ma,np1) = f
    !
    !     Use default program options in DLSEI, and set matrix-vector
    !     printing accuracy parameters.
    !
    prgopt(1) = 1
    idigit = -4
    jdigit = -11
    !
    !     Compute residual norm of known least squares solution.
    !     (to be used to check computed residual norm = RNORML.)
    !
    DO i = 1, ma
      work(i) = DOT_PRODUCT(d(i,1:n),sol) - f(i)
    END DO
    resnrm = NORM2(work(1:ma))
    !
    !     Call DLSEI to get solution in X(*), least squares residual in
    !     RNORML.
    !
    CALL DLSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    !
    !     Compute relative error in problem variable solution and residual
    !     norm computation.
    !
    tnorm = NORM2(sol(1:n))
    err = sol
    CALL DAXPY(n,-1._DP,x,1,err,1)
    cnorm = NORM2(err(1:n))
    relerr = cnorm/tnorm
    relnrm = (resnrm-rnorml)/resnrm
    !
    IF( relerr(1)<=70._DP*SQRT(D1MACH(4)) .AND. relnrm(1)<=5._DP*D1MACH(4) ) THEN
      Ipass = 1
      IF( Kprint>=3 ) WRITE (Lun,99002)
      99002 FORMAT (/' DLSEI PASSED TEST')
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99003) relerr, relnrm
      99003 FORMAT (/' DLSEI FAILED TEST'/' RELERR = ',1P,D20.6/' RELNRM = ',D20.6)
    END IF
    !
    !     Print out known and computed solutions.
    !
    IF( Kprint>=3 ) THEN
      CALL DVOUT(n,err,'('' RESIDUALS FROM KNOWN LEAST SQUARES SOLUTION'')',&
        idigit)
      CALL DVOUT(n,x,'(/'' SOLUTION COMPUTED BY DLSEI'')',jdigit)
    END IF
    !
    IF( Kprint>=2 ) THEN
      IF( Kprint/=2 .OR. Ipass==0 ) THEN
        !
        !           Print out the known and computed residual norms.
        !
        CALL DVOUT(1,resnrm,&
          '(/'' RESIDUAL NORM OF KNOWN LEAST SQUARES SOLUTION'')',jdigit)
        CALL DVOUT(1,rnorml,'(/'' RESIDUAL NORM COMPUTED BY DLSEI'')',jdigit)
        !
        !           Print out the computed solution relative error.
        !
        CALL DVOUT(1,relerr,'(/'' COMPUTED SOLUTION RELATIVE ERROR'')',idigit)
        !
        !           Print out the computed relative error in residual norm.
        !
        CALL DVOUT(1,relnrm,'(/'' COMPUTED RELATIVE ERROR IN RESIDUAL NORM'')'&
          ,idigit)
      END IF
    END IF
    !
    !     Check calls to error processor.
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
    IF( Kprint>=3 ) WRITE (Lun,99004)
    99004 FORMAT (/' 2 ERROR MESSAGES EXPECTED')
    !
    CALL DLSEI(d,0,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    IF( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    END IF
    CALL XERCLR
    !
    prgopt(1) = -1
    CALL DLSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    IF( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    END IF
    CALL XERCLR
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
    CALL XSETF(kontrl)
    IF( fatal ) THEN
      IF( Kprint>=2 ) THEN
        WRITE (Lun,99005)
        99005 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      END IF
    ELSEIF( Kprint>=3 ) THEN
      WRITE (Lun,99006)
      99006 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    END IF
    !
    !     Print PASS/FAIL message.
    !
    IF( Ipass==1 .AND. Kprint>=2 ) WRITE (Lun,99007)
    99007 FORMAT (/' ****************DLSEI PASSED ALL TESTS***************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99008)
    99008 FORMAT (/' ****************DLSEI FAILED SOME TESTS**************')
    RETURN
  END SUBROUTINE DLSEIT
  !** DQCGLS
  SUBROUTINE DQCGLS(Lun,Kprint,Ipass)
    !> Quick check for DGLSS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      DOUBLE PRECISION (QCGLSS-S, DQCGLS-D)
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !      QUICK CHECK SUBROUTINE  DQCGLS  TESTS THE EXECUTION
    !      OF THE GENERAL LINEAR SYSTEM SOLVER, DGLSS .  THE
    !      DGLSS  SUBROUTINE PACKAGE WAS WRITTEN BY T. MANTEUFFEL
    !      (LANL).
    !
    !      A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED
    !      BY DQCGLS.  THE SUMMARY LINE GIVES A COUNT OF THE
    !      NUMBER OF PROBLEMS DETECTED DURING THE TEST.
    !
    !      THE REAL QUANTITIES FOR THE COMPUTED SOLUTION VECTOR
    !      X  AND THE CORRESPONDING  RNORM  ARE COMPARED AGAINST
    !      STORED VALUES.  DISAGREEMENT OCCURS IF A DIFFERENCE
    !      IS SQRT(D1MACH(4) OR MORE.  THE RETURNED VALUE (INTEGER)
    !      OF  INFO  IS ALSO CHECKED.  FOUR CASES ARE RUN, TWO
    !      INVOLVING  LLSIA  AND TWO INVOLVING  ULSIA .
    !
    !      DQCGLS REQUIRES NO INPUT ARGUMENTS.  ON RETURN, NERR
    !      (INTEGER TYPE) CONTAINS THE COUNT OF THE NUMBER OF
    !      PROBLEMS DETECTED BY  QCGLSS .
    !
    !***
    ! **Routines called:**  D1MACH, DGLSS

    !* REVISION HISTORY  (YYMMDD)
    !   811026  DATE WRITTEN
    !   850601  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
    !           including removing an illegal character from column 1, and
    !           editorial changes.  (RWC)
    USE slatec, ONLY : D1MACH, DGLSS
    REAL(DP) :: a(4,4), b(4), delmax, delx, r, rnorm(1), work(50)
    INTEGER :: i, Ipass, j, kk, Kprint, nerr, kprog, kcase, iwork(20), info, Lun
    REAL(DP), PARAMETER :: aa(4,4,2) = RESHAPE( [1._DP, .5_DP, 1._DP, .25_DP, &
      0._DP, 2._DP, 0._DP, 1._DP, 2._DP, -1._DP, 1._DP, 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
      1._DP, 2._DP, -1._DP, 0._DP, 0._DP, 1._DP, 2._DP, 0._DP, -1._DP, 0._DP, 1._DP, 0._DP, &
      1._DP, 0._DP, 1._DP, 0._DP ], [4,4,2] )
    REAL(DP), PARAMETER :: bb(4,2) = RESHAPE( [ 3._DP, 1.5_DP, 2._DP, 1.25_DP, &
      1._DP, 3._DP, 3._DP, 0._DP ], [4,2] )
    REAL(DP), PARAMETER :: xx(4,4) = RESHAPE( [ .9999999999999787_DP, 1.000000000000007_DP, &
      1.000000000000007_DP, 0._DP, .8095238095238102_DP, 1.047619047619044_DP, &
      1.095238095238081_DP, 0._DP, .7777777777777857_DP, 1.444444444444429_DP, &
      .3333333333333393_DP, .5555555555555500_DP, .3333333333333321_DP, 0._DP, &
      -.3333333333333286_DP,.3333333333333286_DP ], [4,4] )
    INTEGER, PARAMETER :: inf(4) = [ 0, 1, 0, 2 ]
    CHARACTER, PARAMETER :: list(2) = [ 'L', 'U' ]
    !* FIRST EXECUTABLE STATEMENT  DQCGLS
    info = 0
    nerr = 0
    r = MAX(SQRT(D1MACH(4)),1.E-12_DP)
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT (/' *  DQCGLS - QUICK CHECK FOR DGLSS (DLLSIA AND DULSIA)'/)
    DO kprog = 1, 2
      DO kcase = 1, 2
        !
        !           FORM BASIC MATRIX  A  AND VECTOR  B .  (CASE 1)
        !
        DO i = 1, 4
          DO j = 1, 4
            a(i,j) = aa(i,j,kprog)
          END DO
          b(i) = bb(i,kprog)
        END DO
        !
        !           MAKE 3 ROWS IDENTICAL FOR CASE 2.
        !
        IF( kcase/=1 ) THEN
          DO i = 2, 3
            DO j = 1, 4
              a(i,j) = a(1,j)
            END DO
            b(i) = b(1)
          END DO
        END IF
        !
        !           SOLVE FOR VECTOR  X .
        !
        info = 0
        IF( kprog==1 ) CALL DGLSS(a,4,4,3,b,4,1,rnorm,work,50,iwork,20,info)
        IF( kprog==2 ) CALL DGLSS(a,4,3,4,b,4,1,rnorm,work,50,iwork,20,info)
        !
        !           TEST COMPUTED  X, RNORM, AND  INFO .
        !
        kk = 2*(kprog-1) + kcase
        delmax = 0._DP
        DO i = 1, 4
          delx = ABS(b(i)-xx(i,kk))
          delmax = MAX(delmax,delx)
        END DO
        !
        IF( Kprint>=3 ) WRITE (Lun,99002) list(kprog), kcase, delmax
        99002 FORMAT (3X,A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',D11.4/)
        IF( delmax>=r ) THEN
          nerr = nerr + 1
          IF( Kprint>=2 ) WRITE (Lun,99003) list(kprog), kcase, delmax
          99003 FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',&
            D11.4/)
        END IF
        !
        IF( Kprint>=3 ) WRITE (Lun,99004) list(kprog), kcase, rnorm
        99004 FORMAT (3X,A,'LSIA, CASE ',I1,'.  RNORM IS ',D11.4/)
        IF( rnorm(1)>=r ) THEN
          nerr = nerr + 1
          IF( Kprint>=2 ) WRITE (Lun,99005) list(kprog), kcase, rnorm
          99005 FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,&
            '.  RNORM (TOO LARGE) IS',D11.4/)
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99006) list(kprog), kcase, info, inf(kk)
        !
        99006 FORMAT (3X,A,'LSIA, CASE ',I1,'.  INFO=',I1,' (SHOULD = ',I1,')'/)
        IF( info/=inf(kk) ) THEN
          nerr = nerr + 1
          IF( Kprint>=2 ) WRITE (Lun,99007) list(kprog), kcase, info, inf(kk)
          99007 FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  INFO=',I1,&
            ' (SHOULD = ',I1,')'/)
        END IF
      END DO
    END DO
    !
    !     SUMMARY PRINT
    !
    Ipass = 0
    IF( nerr==0 ) Ipass = 1
    IF( nerr/=0 .AND. Kprint/=0 ) WRITE (Lun,99008) nerr
    99008 FORMAT (/' **** DQCGLS DETECTED A TOTAL OF ',I2,&
      ' PROBLEMS WITH DGLSS. ****'/)
    IF( nerr==0 .AND. Kprint>1 ) WRITE (Lun,99009)
    99009 FORMAT ('     DQCGLS DETECTED NO PROBLEMS WITH DGLSS.'/)
    RETURN
  END SUBROUTINE DQCGLS
END MODULE TEST28_MOD
!** TEST28
PROGRAM TEST28
  USE TEST28_MOD, ONLY : DLSEIT, DQCGLS
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D5, D9
  !***
  ! **Type:**      DOUBLE PRECISION (TEST28-S, TEST29-D)
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
  !        DLSEI    DGLSS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DLSEIT, DQCGLS, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST28
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
  !     Test DLSEI
  !
  CALL DLSEIT(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test DGLSS
  !
  CALL DQCGLS(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST28 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST28 *************')
  END IF
  STOP
END PROGRAM TEST28
