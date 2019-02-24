MODULE TEST27_MOD
  IMPLICIT NONE

CONTAINS
  !DECK LSEIQX
  SUBROUTINE LSEIQX(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !***BEGIN PROLOGUE  LSEIQX
    !***PURPOSE  Quick check for LSEI.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (LSEIQX-S, DLSEIT-D)
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Hanson, R. J., (SNLA)
    !           Haskell, Karen, (SNLA)
    !***DESCRIPTION
    !
    !   The sample problem solved is from a paper by J. Stoer, in
    !   SIAM Journal of Numerical Analysis, June 1971.
    !
    !***ROUTINES CALLED  LSEI, R1MACH, SAXPY, SCOPY, SDOT, SNRM2, SVOUT
    !***REVISION HISTORY  (YYMMDD)
    !   790216  DATE WRITTEN
    !   890618  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
    !           to use R1MACH(4) rather than R1MACH(3) and cleaned up
    !           FORMATs.  (RWC)
    !   920722  Initialized IP(1) and IP(2) for CALL to LSEI.  (BKS, WRB)
    !   930214  Declarations sections added, code revised to test error
    !           returns for all values of KPRINT and code polished.  (WRB)
    !***END PROLOGUE  LSEIQX
    !     .. Scalar Arguments ..
    INTEGER Ipass, Kprint, Lun
    !     .. Local Scalars ..
    REAL cnorm, relerr, relnrm, resnrm, rnorme, rnorml, tnorm
    INTEGER i, idigit, jdigit, kontrl, ma, mdd, me, meap1, mep1, mg, &
      mode, n, nerr, np1
    LOGICAL fatal
    !     .. Local Arrays ..
    REAL a(6,5), d(11,6), err(5), f(6), g(5,5), h(5), prgopt(4), sol(5)&
      , work(105), x(5)
    INTEGER ip(17)
    !     .. External Functions ..
    REAL R1MACH, SDOT, SNRM2
    INTEGER NUMXER
    EXTERNAL NUMXER, R1MACH, SDOT, SNRM2
    !     .. External Subroutines ..
    EXTERNAL LSEI, SAXPY, SCOPY, SVOUT, XGETF, XSETF
    !     .. Intrinsic Functions ..
    INTRINSIC SQRT
    !     .. Data statements ..
    !
    !     Define the data arrays for the example.  The array A contains
    !     the least squares equations.  (There are no equality constraints
    !     in this example).
    !
    DATA a(1,1), a(1,2), a(1,3), a(1,4), a(1,5)/ - 74., 80., 18., &
      -11., -4./
    DATA a(2,1), a(2,2), a(2,3), a(2,4), a(2,5)/14., -69., 21., 28., &
      0./
    DATA a(3,1), a(3,2), a(3,3), a(3,4), a(3,5)/66., -72., -5., 7., &
      1./
    DATA a(4,1), a(4,2), a(4,3), a(4,4), a(4,5)/ - 12., 66., -30., &
      -23., 3./
    DATA a(5,1), a(5,2), a(5,3), a(5,4), a(5,5)/3., 8., -7., -4., 1./
    DATA a(6,1), a(6,2), a(6,3), a(6,4), a(6,5)/4., -12., 4., 4., 0./
    !
    !     The array G contains the inequality constraint equations,
    !     written in the sense
    !     (row vector)*(solution vector) .GE. (given value).
    !
    DATA g(1,1), g(1,2), g(1,3), g(1,4), g(1,5)/ - 1., -1., -1., -1., &
      -1./
    DATA g(2,1), g(2,2), g(2,3), g(2,4), g(2,5)/10., 10., -3., 5., 4./
    DATA g(3,1), g(3,2), g(3,3), g(3,4), g(3,5)/ - 8., 1., -2., -5., &
      3./
    DATA g(4,1), g(4,2), g(4,3), g(4,4), g(4,5)/8., -1., 2., 5., -3./
    DATA g(5,1), g(5,2), g(5,3), g(5,4), g(5,5)/ - 4., -2., 3., -5., &
      1./
    !
    !     Define the least squares right-side vector.
    !
    DATA f(1), f(2), f(3), f(4), f(5), f(6)/ - 5., -9., 708., 4165., &
      -13266., 8409./
    !
    !     Define the inequality constraint right-side vector.
    !
    DATA h(1), h(2), h(3), h(4), h(5)/ - 5., 20., -40., 11., -30./
    !
    !     Define the vector that is the known solution.
    !
    DATA sol(1), sol(2), sol(3), sol(4), sol(5)/1., 2., -1., 3., -4./
    !***FIRST EXECUTABLE STATEMENT  LSEIQX
    IF ( Kprint>=2 ) WRITE (Lun,99001)
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
    DO i = 1, n
      !
      !        Copy the i-th column of the inequality constraint matrix into
      !        the work array.
      !
      CALL SCOPY(mg,g(1,i),1,d(meap1,i),1)
      !
      !        Copy the i-th column of the least squares matrix into the work
      !        array.
      !
      CALL SCOPY(ma,a(1,i),1,d(mep1,i),1)
    ENDDO
    !
    !     Copy the right-side vectors into the work array in compatible
    !     order.
    !
    CALL SCOPY(mg,h,1,d(meap1,np1),1)
    CALL SCOPY(ma,f,1,d(mep1,np1),1)
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
      work(i) = SDOT(n,d(i,1),mdd,sol,1) - f(i)
    ENDDO
    resnrm = SNRM2(ma,work,1)
    !
    !     Call LSEI to get solution in X(*), least squares residual in
    !     RNORML.
    !
    CALL LSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
    !
    !     Compute relative error in problem variable solution and residual
    !     norm computation.
    !
    tnorm = SNRM2(n,sol,1)
    CALL SCOPY(n,sol,1,err,1)
    CALL SAXPY(n,-1.0E0,x,1,err,1)
    cnorm = SNRM2(n,err,1)
    relerr = cnorm/tnorm
    relnrm = (resnrm-rnorml)/resnrm
    !
    IF ( relerr<=70.0E0*SQRT(R1MACH(4)).AND.relnrm<=5.0E0*R1MACH(4) ) THEN
      Ipass = 1
      IF ( Kprint>=3 ) WRITE (Lun,99002)
      99002   FORMAT (/' LSEI PASSED TEST')
    ELSE
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (Lun,99003) relerr, relnrm
      99003   FORMAT (/' LSEI FAILED TEST'/' RELERR = ',1P,E20.6/' RELNRM = ',E20.6)
    ENDIF
    !
    !     Print out known and computed solutions.
    !
    IF ( Kprint>=3 ) THEN
      CALL SVOUT(n,err,'('' RESIDUALS FROM KNOWN LEAST SQUARES SOLUTION'')',&
        idigit)
      CALL SVOUT(n,x,'(/'' SOLUTION COMPUTED BY LSEI'')',jdigit)
    ENDIF
    !
    IF ( Kprint>=2 ) THEN
      IF ( Kprint/=2.OR.Ipass==0 ) THEN
        !
        !           Print out the known and computed residual norms.
        !
        CALL SVOUT(1,resnrm,&
          '(/'' RESIDUAL NORM OF KNOWN LEAST SQUARES SOLUTION'')',&
          jdigit)
        CALL SVOUT(1,rnorml,'(/'' RESIDUAL NORM COMPUTED BY LSEI'')',jdigit)
        !
        !           Print out the computed solution relative error.
        !
        CALL SVOUT(1,relerr,'(/'' COMPUTED SOLUTION RELATIVE ERROR'')',idigit)
        !
        !           Print out the computed relative error in residual norm.
        !
        CALL SVOUT(1,relnrm,'(/'' COMPUTED RELATIVE ERROR IN RESIDUAL NORM'')'&
          ,idigit)
      ENDIF
    ENDIF
    !
    !     Check calls to error processor.
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
    IF ( Kprint>=3 ) WRITE (Lun,99004)
    99004 FORMAT (/' 2 ERROR MESSAGES EXPECTED')
    !
    CALL LSEI(d,0,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
    IF ( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    !
    prgopt(1) = -1
    CALL LSEI(d,mdd,me,ma,mg,n,prgopt,x,rnorme,rnorml,mode,work,ip)
    IF ( NUMXER(nerr)/=2 ) THEN
      Ipass = 0
      fatal = .TRUE.
    ENDIF
    CALL XERCLR
    !
    !     Restore KONTRL and check to see if the tests of error detection
    !     passed.
    !
    CALL XSETF(kontrl)
    IF ( fatal ) THEN
      IF ( Kprint>=2 ) THEN
        WRITE (Lun,99005)
        99005     FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
      ENDIF
    ELSEIF ( Kprint>=3 ) THEN
      WRITE (Lun,99006)
      99006   FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
    ENDIF
    !
    !     Print PASS/FAIL message.
    !
    IF ( Ipass==1.AND.Kprint>=2 ) WRITE (Lun,99007)
    99007 FORMAT (/' ****************LSEI PASSED ALL TESTS***************')
    IF ( Ipass==0.AND.Kprint>=1 ) WRITE (Lun,99008)
    99008 FORMAT (/' ****************LSEI FAILED SOME TESTS**************')
    RETURN
  END SUBROUTINE LSEIQX
  !DECK QCGLSS
  SUBROUTINE QCGLSS(Lun,Kprint,Ipass)
    IMPLICIT NONE
    INTEGER i, Ipass, j, kk, Kprint
    REAL R1MACH, rnorm
    !***BEGIN PROLOGUE  QCGLSS
    !***PURPOSE  Quick check for SGLSS.
    !***LIBRARY   SLATEC
    !***TYPE      SINGLE PRECISION (QCGLSS-S, DQCGLS-D)
    !***AUTHOR  Voorhees, E. A., (LANL)
    !***DESCRIPTION
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
    !      IS SQRT(R1MACH(4) OR MORE.  THE RETURNED VALUE (INTEGER)
    !      OF  INFO  IS ALSO CHECKED.  FOUR CASES ARE RUN, TWO
    !      INVOLVING  LLSIA  AND TWO INVOLVING  ULSIA .
    !
    !      QCGLSS REQUIRES NO INPUT ARGUMENTS.  ON RETURN, NERR
    !      (INTEGER TYPE) CONTAINS THE COUNT OF THE NUMBER OF
    !      PROBLEMS DETECTED BY  QCGLSS .
    !
    !***ROUTINES CALLED  R1MACH, SGLSS
    !***REVISION HISTORY  (YYMMDD)
    !   811026  DATE WRITTEN
    !   820801  REVISION DATE from Version 3.2
    !   891214  Prologue converted to Version 4.0 format.  (BAB)
    !   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs,
    !           including removing an illegal character from column 1, and
    !           editorial changes.  (RWC)
    !***END PROLOGUE  QCGLSS
    REAL aa(4,4,2), a(4,4), bb(4,2), b(4), xx(4,4), delmax, delx, r
    REAL work(20)
    CHARACTER :: list(2)
    INTEGER inf(4), nerr, kprog, kcase
    INTEGER iwork(7), info, Lun
    DATA aa/1., .5, 1., .25, 0., 2., 0., 1., 2., -1., 1., 0., 0. ,&
      0., 0., 0., 1., 2., -1., 0., 0., 1., 2., 0., -1., 0. ,&
      1., 0., 1., 0., 1., 0./
    DATA bb/3., 1.5, 2., 1.25, 1., 3., 3., 0./
    DATA xx/.9999999999999787, 1.000000000000007, 1.000000000000007, 0. ,&
      .8095238095238102, 1.047619047619044, 1.095238095238081, 0. ,&
      .7777777777777857, 1.444444444444429, .3333333333333393 ,&
      .5555555555555500, .3333333333333321, 0.0, -.3333333333333286 ,&
      .3333333333333286/
    DATA inf/0, 1, 0, 2/
    DATA list/'L', 'U'/
    !***FIRST EXECUTABLE STATEMENT  QCGLSS
    info = 0
    nerr = 0
    r = SQRT(R1MACH(4))
    IF ( Kprint>=2 ) WRITE (Lun,99001)
    99001 FORMAT (/' *    QCGLSS - QUICK CHECK FOR SGLSS (LLSIA AND ULSIA)'/)
    DO kprog = 1, 2
      DO kcase = 1, 2
        !
        !           FORM BASIC MATRIX  A  AND VECTOR  B .  (CASE 1)
        !
        DO i = 1, 4
          DO j = 1, 4
            a(i,j) = aa(i,j,kprog)
          ENDDO
          b(i) = bb(i,kprog)
        ENDDO
        !
        !           MAKE 3 ROWS IDENTICAL FOR CASE 2.
        !
        IF ( kcase/=1 ) THEN
          DO i = 2, 3
            DO j = 1, 4
              a(i,j) = a(1,j)
            ENDDO
            b(i) = b(1)
          ENDDO
        ENDIF
        !
        !           SOLVE FOR VECTOR  X .
        !
        info = 0
        IF ( kprog==1 ) CALL SGLSS(a,4,4,3,b,4,1,rnorm,work,20,iwork,7,info)
        IF ( kprog==2 ) CALL SGLSS(a,4,3,4,b,4,1,rnorm,work,20,iwork,7,info)
        !
        !           TEST COMPUTED  X, RNORM, AND  INFO .
        !
        kk = 2*(kprog-1) + kcase
        delmax = 0.0E0
        DO i = 1, 4
          delx = ABS(b(i)-xx(i,kk))
          delmax = MAX(delmax,delx)
        ENDDO
        !
        IF ( Kprint>=3 ) WRITE (Lun,99002) list(kprog), kcase, delmax
        !
        99002     FORMAT (3X,A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',E11.4/)
        IF ( delmax>=r ) THEN
          nerr = nerr + 1
          IF ( Kprint>=2 ) WRITE (Lun,99003) list(kprog), kcase, delmax
          99003       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  MAX ABS ERROR OF',&
            E11.4/)
        ENDIF
        IF ( Kprint>=3 ) WRITE (Lun,99004) list(kprog), kcase, rnorm
        99004     FORMAT (3X,A,'LSIA, CASE ',I1,'.  RNORM IS ',E11.4/)
        IF ( rnorm>r ) THEN
          nerr = nerr + 1
          IF ( Kprint>=2 ) WRITE (Lun,99005) list(kprog), kcase, rnorm
          99005       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,&
            '.  RNORM (TOO LARGE) IS',E11.4/)
        ENDIF
        !
        IF ( Kprint>=3 ) WRITE (Lun,99006) list(kprog), kcase, info ,&
          inf(kk)
        99006     FORMAT (3X,A,'LSIA, CASE ',I1,'.  INFO=',I1,' (SHOULD = ',I1,')'/)
        IF ( info/=inf(kk) ) THEN
          nerr = nerr + 1
          IF ( Kprint>=2 ) WRITE (Lun,99007) list(kprog), kcase, info ,&
            inf(kk)
          99007       FORMAT ('   PROBLEM WITH ',A,'LSIA, CASE ',I1,'.  INFO=',I1,&
            ' (SHOULD = ',I1,')'/)
        ENDIF
      ENDDO
    ENDDO
    !
    !     SUMMARY PRINT
    !
    Ipass = 0
    IF ( nerr==0 ) Ipass = 1
    IF ( nerr/=0.AND.Kprint/=0 ) WRITE (Lun,99008) nerr
    99008 FORMAT (/' **** QCGLSS DETECTED A TOTAL OF ',I2,&
      ' PROBLEMS WITH SGLSS. ****'/)
    IF ( nerr==0.AND.Kprint>1 ) WRITE (Lun,99009)
    99009 FORMAT ('     QCGLSS DETECTED NO PROBLEMS WITH SGLSS.'/)
    RETURN
  END SUBROUTINE QCGLSS
END MODULE TEST27_MOD
!DECK TEST27
PROGRAM TEST27
  USE TEST27_MOD
  IMPLICIT NONE
  INTEGER I1MACH
  !***BEGIN PROLOGUE  TEST27
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  D5, D9
  !***TYPE      SINGLE PRECISION (TEST27-S, TEST28-D)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  ! *Arguments:
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
  ! *Description:
  !     Driver for testing SLATEC subprograms
  !        LSEI     SGLSS
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, LSEIQX, QCGLSS, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST27
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST27
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (lin,'(I1)') kprint
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test LSEI
  !
  CALL LSEIQX(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test SGLSS
  !
  CALL QCGLSS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001   FORMAT (/' --------------TEST27 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002   FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST27 *************')
  ENDIF
  STOP
END PROGRAM TEST27
