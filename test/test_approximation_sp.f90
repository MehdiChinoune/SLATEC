MODULE TEST27_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** LSEIQX
  SUBROUTINE LSEIQX(Lun,Kprint,Ipass)
    !> Quick check for LSEI.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (LSEIQX-S, DLSEIT-D)
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
    ! **Routines called:**  LSEI, R1MACH, SAXPY, SCOPY, SDOT, SNRM2, SVOUT

    !* REVISION HISTORY  (YYMMDD)
    !   790216  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
    !           to use eps_sp rather than eps_2_sp and cleaned up FORMATs.  (RWC)
    !   920722  Initialized IP(1) and IP(2) for CALL to LSEI.  (BKS, WRB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    USE service, ONLY : eps_sp, SVOUT
    USE approximation, ONLY : LSEI
    USE blas, ONLY :  SAXPY
    !     .. Scalar Arguments ..
    INTEGER :: Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL(SP) :: cnorm, relerr(1), relnrm(1), resnrm(1), rnorme, rnorml(1), tnorm
    INTEGER :: i, idigit, jdigit, ma, mdd, me, meap1, mep1, mg, &
      mode, n, np1
    LOGICAL :: fatal
    !     .. Local Arrays ..
    REAL(SP) :: d(11,6), err(5), prgopt(4), work(105), x(5)
    INTEGER :: ip(17)
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !     .. Data statements ..
    !
    !     Define the data arrays for the example.  The array A contains
    !     the least squares equations.  (There are no equality constraints
    !     in this example).
    !
    REAL(SP), PARAMETER :: a(6,5) = RESHAPE( [ &
      -74._SP, 80._SP, 18._SP, -11._SP, -4._SP,    14._SP, -69._SP, 21._SP, 28._SP, 0._SP, &
      66._SP, -72._SP, -5._SP, 7._SP, 1._SP,      -12._SP, 66._SP, -30._SP, -23._SP, 3._SP, &
      3._SP, 8._SP, -7._SP, -4._SP, 1._SP,         4._SP, -12._SP, 4._SP, 4._SP, 0._SP ], &
      [6,5], ORDER = [2,1] )
    !
    !     The array G contains the inequality constraint equations,
    !     written in the sense
    !     (row vector)*(solution vector) >= (given value).
    !
    REAL(SP), PARAMETER :: g(5,5) = RESHAPE( [ -1._SP, -1._SP, -1._SP, -1._SP, -1._SP, &
      10._SP, 10._SP, -3._SP, 5._SP, 4._SP,    -8._SP, 1._SP, -2._SP, -5._SP, 3._SP, &
      8._SP, -1._SP, 2._SP, 5._SP, -3._SP,     -4._SP, -2._SP, 3._SP, -5._SP, 1._SP ], &
      [5,5], ORDER = [2,1] )
    !
    !     Define the least squares right-side vector.
    !
    REAL(SP), PARAMETER :: f(6) = [ -5._SP, -9._SP, 708._SP, 4165._SP, &
      -13266._SP, 8409._SP ]
    !
    !     Define the inequality constraint right-side vector.
    !
    REAL(SP), PARAMETER :: h(5) = [ -5._SP, 20._SP, -40._SP, 11._SP, -30._SP ]
    !
    !     Define the vector that is the known solution.
    !
    REAL(SP), PARAMETER :: sol(5) = [ 1._SP, 2._SP, -1._SP, 3._SP, -4._SP ]
    !* FIRST EXECUTABLE STATEMENT  LSEIQX
    IF( Kprint>=2 ) WRITE (Lun,99001)
    !
    99001 FORMAT ('1TEST OF SUBROUTINE LSEI')
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
    !     Copy the right-side vectors into the work array in compatible
    !     order.
    !
    d(ma+1:mdd,np1) = h
    d(1:ma,np1) = f
    !
    !     Use default program options in LSEI, and set matrix-vector
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
    !     Call LSEI to get solution in X(*), least squares residual in
    !     RNORML.
    !
    CALL LSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
    !
    !     Compute relative error in problem variable solution and residual
    !     norm computation.
    !
    tnorm = NORM2(sol(1:n))
    err = sol
    CALL SAXPY(n,-1._SP,x,1,err,1)
    cnorm = NORM2(err(1:n))
    relerr = cnorm/tnorm
    relnrm = (resnrm-rnorml)/resnrm
    !
    IF( relerr(1)<=70._SP*SQRT(eps_sp) .AND. relnrm(1)<=5._SP*eps_sp ) THEN
      Ipass = 1
      IF( Kprint>=3 ) WRITE (Lun,99002)
      99002 FORMAT (/' LSEI PASSED TEST')
    ELSE
      Ipass = 0
      IF( Kprint>=2 ) WRITE (Lun,99003) relerr, relnrm
      99003 FORMAT (/' LSEI FAILED TEST'/' RELERR = ',1P,E20.6/' RELNRM = ',E20.6)
    END IF
    !
    !     Print out known and computed solutions.
    !
    IF( Kprint>=3 ) THEN
      CALL SVOUT(n,err,'('' RESIDUALS FROM KNOWN LEAST SQUARES SOLUTION'')',idigit)
      CALL SVOUT(n,x,'(/'' SOLUTION COMPUTED BY LSEI'')',jdigit)
    END IF
    !
    IF( Kprint>=2 ) THEN
      IF( Kprint/=2 .OR. Ipass==0 ) THEN
        !
        !           Print out the known and computed residual norms.
        !
        CALL SVOUT(1,resnrm,'(/'' RESIDUAL NORM OF KNOWN LEAST SQUARES SOLUTION'')',jdigit)
        CALL SVOUT(1,rnorml,'(/'' RESIDUAL NORM COMPUTED BY LSEI'')',jdigit)
        !
        !           Print out the computed solution relative error.
        !
        CALL SVOUT(1,relerr,'(/'' COMPUTED SOLUTION RELATIVE ERROR'')',idigit)
        !
        !           Print out the computed relative error in residual norm.
        !
        CALL SVOUT(1,relnrm,'(/'' COMPUTED RELATIVE ERROR IN RESIDUAL NORM'')',idigit)
      END IF
    END IF
    !
    !     Check calls to error processor.
    !
    fatal = .FALSE.
!    IF( Kprint<=2 ) THEN
!      control_xer = 0
!    ELSE
!      control_xer = 1
!    END IF
!    num_xer = 0
    !
!    IF( Kprint>=3 ) WRITE (Lun,99004)
!    99004 FORMAT (/' 2 ERROR MESSAGES EXPECTED')
    !
!    CALL LSEI(d,0,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
    !
!    prgopt(1) = -1
!    CALL LSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml(1),mode,work,ip)
!    IF( num_xer/=2 ) THEN
!      Ipass = 0
!      fatal = .TRUE.
!    END IF
!    num_xer = 0
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
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
    99007 FORMAT (/' ****************LSEI PASSED ALL TESTS***************')
    IF( Ipass==0 .AND. Kprint>=1 ) WRITE (Lun,99008)
    99008 FORMAT (/' ****************LSEI FAILED SOME TESTS**************')
    RETURN
  END SUBROUTINE LSEIQX
  !** QCGLSS
  SUBROUTINE QCGLSS(Lun,Kprint,Ipass)
    !> Quick check for SGLSS.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Type:**      SINGLE PRECISION (QCGLSS-S, DQCGLS-D)
    !***
    ! **Author:**  Voorhees, E. A., (LANL)
    !***
    ! **Description:**
    !
    !      QUICK CHECK SUBROUTINE  QCGLSS  TESTS THE EXECUTION
    !      OF THE GENERAL LINEAR SYSTEM SOLVER, SGLSS .  THE
    !      SGLSS  SUBROUTINE PACKAGE WAS WRITTEN BY T. MANTEUFFEL
    !      (LANL).
    !
    !      A TITLE LINE AND A SUMMARY LINE ARE ALWAYS OUTPUTTED
    !      BY QCGLSS.  THE SUMMARY LINE GIVES A COUNT OF THE
    !      NUMBER OF PROBLEMS DETECTED DURING THE TEST.
    !
    !      THE REAL QUANTITIES FOR THE COMPUTED SOLUTION VECTOR
    !      X  AND THE CORRESPONDING  RNORM  ARE COMPARED AGAINST
    !      STORED VALUES.  DISAGREEMENT OCCURS IF A DIFFERENCE
    !      IS SQRT(eps_sp OR MORE.  THE RETURNED VALUE (INTEGER)
    !      OF  INFO  IS ALSO CHECKED.  FOUR CASES ARE RUN, TWO
    !      INVOLVING  LLSIA  AND TWO INVOLVING  ULSIA .
    !
    !      QCGLSS REQUIRES NO INPUT ARGUMENTS.  ON RETURN, NERR
    !      (INTEGER TYPE) CONTAINS THE COUNT OF THE NUMBER OF
    !      PROBLEMS DETECTED BY  QCGLSS .
    !
    !***
    ! **Routines called:**  R1MACH, SGLSS

    !* REVISION HISTORY  (YYMMDD)
    !   811026  DATE WRITTEN
    !   820801  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs, including
    !           removing an illegal character from column 1, and editorial changes.  (RWC)
    USE service, ONLY : eps_sp
    USE linear, ONLY : SGLSS
    !
    INTEGER :: i, Ipass, j, kk, Kprint, nerr, kprog, kcase, iwork(7), info, Lun
    REAL(SP) :: rnorm(1), a(4,4), b(4), delmax, delx, r, work(20)
    REAL(SP), PARAMETER :: aa(4,4,2) = RESHAPE( [ 1._SP, .5_SP, 1._SP, .25_SP, &
      0._SP, 2._SP, 0._SP, 1._SP,   2._SP, -1._SP, 1._SP, 0._SP, &
      0._SP, 0._SP, 0._SP, 0._SP,   1._SP, 2._SP, -1._SP, 0._SP, &
      0._SP, 1._SP, 2._SP, 0._SP,  -1._SP, 0._SP, 1._SP, 0._SP, &
      1._SP, 0._SP, 1._SP, 0._SP ], [4,4,2] )
    REAL(SP), PARAMETER :: bb(4,2) = RESHAPE( [ 3._SP, 1.5_SP, 2._SP, 1.25_SP, &
      1._SP, 3._SP, 3._SP, 0._SP ], [4,2] )
    REAL(SP), PARAMETER :: xx(4,4) = RESHAPE( [ &
      .9999999999999787_SP, 1.000000000000007_SP, 1.000000000000007_SP, 0._SP, &
      .8095238095238102_SP, 1.047619047619044_SP, 1.095238095238081_SP, 0._SP, &
      .7777777777777857_SP, 1.444444444444429_SP, .3333333333333393_SP, .5555555555555500_SP, &
      .3333333333333321_SP, 0.0_SP, -.3333333333333286_SP, .3333333333333286_SP ], [4,4] )
    INTEGER, PARAMETER :: inf(4) = [ 0, 1, 0, 2 ]
    CHARACTER, PARAMETER :: list(2) = [ 'L', 'U' ]
    !* FIRST EXECUTABLE STATEMENT  QCGLSS
    info = 0
    nerr = 0
    r = SQRT(eps_sp)
    IF( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT (/' *    QCGLSS - QUICK CHECK FOR SGLSS (LLSIA AND ULSIA)'/)
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
        IF( kprog==1 ) CALL SGLSS(a,4,4,3,b,4,1,rnorm,work,20,iwork,7,info)
        IF( kprog==2 ) CALL SGLSS(a,4,3,4,b,4,1,rnorm,work,20,iwork,7,info)
        !
        !           TEST COMPUTED  X, RNORM, AND  INFO .
        !
        kk = 2*(kprog-1) + kcase
        delmax = 0._SP
        DO i = 1, 4
          delx = ABS(b(i)-xx(i,kk))
          delmax = MAX(delmax,delx)
        END DO
        !
        IF( Kprint>=3 ) WRITE (Lun,99002) list(kprog), kcase, delmax
        !
        99002 FORMAT (3X,A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',E11.4/)
        IF( delmax>=r ) THEN
          nerr = nerr + 1
          IF( Kprint>=2 ) WRITE (Lun,99003) list(kprog), kcase, delmax
          99003 FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',&
            E11.4/)
        END IF
        IF( Kprint>=3 ) WRITE (Lun,99004) list(kprog), kcase, rnorm
        99004 FORMAT (3X,A,'LSIA, CASE ',I1,'.  RNORM IS ',E11.4/)
        IF( rnorm(1)>r ) THEN
          nerr = nerr + 1
          IF( Kprint>=2 ) WRITE (Lun,99005) list(kprog), kcase, rnorm
          99005 FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,&
            '.  RNORM (TOO LARGE) IS',E11.4/)
        END IF
        !
        IF( Kprint>=3 ) WRITE (Lun,99006) list(kprog), kcase, info, inf(kk)
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
    99008 FORMAT (/' **** QCGLSS DETECTED A TOTAL OF ',I2,&
      ' PROBLEMS WITH SGLSS. ****'/)
    IF( nerr==0 .AND. Kprint>1 ) WRITE (Lun,99009)
    99009 FORMAT ('     QCGLSS DETECTED NO PROBLEMS WITH SGLSS.'/)
    RETURN
  END SUBROUTINE QCGLSS
END MODULE TEST27_MOD
!** TEST27
PROGRAM TEST27
  USE TEST27_MOD, ONLY : LSEIQX, QCGLSS
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D5, D9
  !***
  ! **Type:**      SINGLE PRECISION (TEST27-S, TEST28-D)
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
  !        LSEI     SGLSS
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, LSEIQX, QCGLSS, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST27
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test LSEI
  !
  CALL LSEIQX(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Test SGLSS
  !
  CALL QCGLSS(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST27 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST27 *************')
  END IF
  STOP
END PROGRAM TEST27
