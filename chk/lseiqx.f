*DECK LSEIQX
      SUBROUTINE LSEIQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  LSEIQX
C***PURPOSE  Quick check for LSEI.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSEIQX-S, DLSEIT-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Hanson, R. J., (SNLA)
C           Haskell, Karen, (SNLA)
C***DESCRIPTION
C
C   The sample problem solved is from a paper by J. Stoer, in
C   SIAM Journal of Numerical Analysis, June 1971.
C
C***ROUTINES CALLED  LSEI, R1MACH, SAXPY, SCOPY, SDOT, SNRM2, SVOUT
C***REVISION HISTORY  (YYMMDD)
C   790216  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Restructured using IF-THEN-ELSE-ENDIF, modified tolerances
C           to use R1MACH(4) rather than R1MACH(3) and cleaned up
C           FORMATs.  (RWC)
C   920722  Initialized IP(1) and IP(2) for CALL to LSEI.  (BKS, WRB)
C   930214  Declarations sections added, code revised to test error
C           returns for all values of KPRINT and code polished.  (WRB)
C***END PROLOGUE  LSEIQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL CNORM, RELERR, RELNRM, RESNRM, RNORME, RNORML, TNORM
      INTEGER I, IDIGIT, JDIGIT, KONTRL, MA, MDD, ME, MEAP1, MEP1, MG,
     *        MODE, N, NERR, NP1
      LOGICAL FATAL
C     .. Local Arrays ..
      REAL A(6,5), D(11,6), ERR(5), F(6), G(5,5), H(5), PRGOPT(4),
     *     SOL(5), WORK(105), X(5)
      INTEGER IP(17)
C     .. External Functions ..
      REAL R1MACH, SDOT, SNRM2
      INTEGER NUMXER
      EXTERNAL NUMXER, R1MACH, SDOT, SNRM2
C     .. External Subroutines ..
      EXTERNAL LSEI, SAXPY, SCOPY, SVOUT, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     .. Data statements ..
C
C     Define the data arrays for the example.  The array A contains
C     the least squares equations.  (There are no equality constraints
C     in this example).
C
      DATA A(1,1),A(1,2),A(1,3),A(1,4),A(1,5)
     *     /-74.,80.,18.,-11.,-4./
      DATA A(2,1),A(2,2),A(2,3),A(2,4),A(2,5)
     *     /14.,-69.,21.,28.,0./
      DATA A(3,1),A(3,2),A(3,3),A(3,4),A(3,5)
     *     /66.,-72.,-5.,7.,1./
      DATA A(4,1),A(4,2),A(4,3),A(4,4),A(4,5)
     *     /-12.,66.,-30.,-23.,3./
      DATA A(5,1),A(5,2),A(5,3),A(5,4),A(5,5)
     *     /3.,8.,-7.,-4.,1./
      DATA A(6,1),A(6,2),A(6,3),A(6,4),A(6,5)
     *     /4.,-12.,4.,4.,0./
C
C     The array G contains the inequality constraint equations,
C     written in the sense
C     (row vector)*(solution vector) .GE. (given value).
C
      DATA G(1,1),G(1,2),G(1,3),G(1,4),G(1,5)
     *     /-1.,-1.,-1.,-1.,-1./
      DATA G(2,1),G(2,2),G(2,3),G(2,4),G(2,5)
     *     /10.,10.,-3.,5.,4./
      DATA G(3,1),G(3,2),G(3,3),G(3,4),G(3,5)
     *     /-8.,1.,-2.,-5.,3./
      DATA G(4,1),G(4,2),G(4,3),G(4,4),G(4,5)
     *     /8.,-1.,2.,5.,-3./
      DATA G(5,1),G(5,2),G(5,3),G(5,4),G(5,5)
     *     /-4.,-2.,3.,-5.,1./
C
C     Define the least squares right-side vector.
C
      DATA F(1),F(2),F(3),F(4),F(5),F(6)
     *     /-5.,-9.,708.,4165.,-13266.,8409./
C
C     Define the inequality constraint right-side vector.
C
      DATA H(1),H(2),H(3),H(4),H(5)
     *     /-5.,20.,-40.,11.,-30./
C
C     Define the vector that is the known solution.
C
      DATA SOL(1),SOL(2),SOL(3),SOL(4),SOL(5)
     *     /1.,2.,-1.,3.,-4./
C***FIRST EXECUTABLE STATEMENT  LSEIQX
      IF (KPRINT .GE. 2) WRITE (LUN, 9000)
C
C     Define the matrix dimensions, number of least squares equations,
C     number of equality constraints, total number of equations, and
C     number of variables.  Set ME=0 to indicate there are no equality
C     constraints.
C
      MDD = 11
      MA = 6
      MG = 5
      N = 5
      ME = 0
C
      IP(1) = 105
      IP(2) = 17
C
      NP1 = N + 1
      MEP1 = ME + 1
      MEAP1 = ME + MA + 1
C
C     Copy the problem matrices.
C
      DO 10 I = 1, N
C
C        Copy the i-th column of the inequality constraint matrix into
C        the work array.
C
         CALL SCOPY(MG, G(1,I), 1, D(MEAP1,I), 1)
C
C        Copy the i-th column of the least squares matrix into the work
C        array.
C
         CALL SCOPY(MA, A(1,I), 1, D(MEP1,I), 1)
   10 CONTINUE
C
C     Copy the right-side vectors into the work array in compatible
C     order.
C
      CALL SCOPY(MG, H, 1, D(MEAP1,NP1), 1)
      CALL SCOPY(MA, F, 1, D(MEP1,NP1),  1)
C
C     Use default program options in LSEI, and set matrix-vector
C     printing accuracy parameters.
C
      PRGOPT(1) = 1
      IDIGIT = -4
      JDIGIT = -11
C
C     Compute residual norm of known least squares solution.
C     (to be used to check computed residual norm = RNORML.)
C
      DO 20 I = 1, MA
         WORK(I) = SDOT(N,D(I,1),MDD,SOL,1) - F(I)
   20 CONTINUE
      RESNRM = SNRM2(MA,WORK,1)
C
C     Call LSEI to get solution in X(*), least squares residual in
C     RNORML.
C
      CALL LSEI(D, MDD, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML, MODE,
     *   WORK, IP)
C
C     Compute relative error in problem variable solution and residual
C     norm computation.
C
      TNORM = SNRM2(N,SOL,1)
      CALL SCOPY(N, SOL, 1, ERR, 1)
      CALL SAXPY(N, -1.0E0, X, 1, ERR, 1)
      CNORM = SNRM2(N, ERR, 1)
      RELERR = CNORM/TNORM
      RELNRM = (RESNRM-RNORML)/RESNRM
C
      IF (RELERR .LE. 70.0E0*SQRT(R1MACH(4)) .AND.
     *    RELNRM .LE.  5.0E0*R1MACH(4)) THEN
         IPASS = 1
         IF (KPRINT .GE. 3) WRITE (LUN, 9010)
      ELSE
         IPASS = 0
         IF (KPRINT .GE. 2) WRITE (LUN, 9020) RELERR, RELNRM
      ENDIF
C
C     Print out known and computed solutions.
C
      IF (KPRINT .GE. 3) THEN
         CALL SVOUT(N, ERR,
     *      '('' RESIDUALS FROM KNOWN LEAST SQUARES SOLUTION'')',
     *      IDIGIT)
         CALL SVOUT(N, X, '(/'' SOLUTION COMPUTED BY LSEI'')', JDIGIT)
      ENDIF
C
      IF (KPRINT .GE. 2) THEN
         IF (.NOT.(KPRINT.EQ.2 .AND. IPASS.NE.0)) THEN
C
C           Print out the known and computed residual norms.
C
            CALL SVOUT(1, RESNRM,
     *         '(/'' RESIDUAL NORM OF KNOWN LEAST SQUARES SOLUTION'')',
     *         JDIGIT)
            CALL SVOUT(1, RNORML,
     *         '(/'' RESIDUAL NORM COMPUTED BY LSEI'')', JDIGIT)
C
C           Print out the computed solution relative error.
C
            CALL SVOUT(1, RELERR,
     *         '(/'' COMPUTED SOLUTION RELATIVE ERROR'')', IDIGIT)
C
C           Print out the computed relative error in residual norm.
C
            CALL SVOUT(1, RELNRM,
     *       '(/'' COMPUTED RELATIVE ERROR IN RESIDUAL NORM'')', IDIGIT)
         ENDIF
      ENDIF
C
C     Check calls to error processor.
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
      FATAL = .FALSE.
      CALL XERCLR
C
      IF (KPRINT .GE. 3) WRITE (LUN, 9030)
C
      CALL LSEI (D, 0, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML,
     *           MODE, WORK, IP)
      IF (NUMXER(NERR) .NE. 2) THEN
         IPASS = 0
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      PRGOPT(1) = -1
      CALL LSEI (D, MDD, ME, MA, MG, N, PRGOPT, X, RNORME, RNORML,
     *           MODE, WORK, IP)
      IF (NUMXER(NERR) .NE. 2) THEN
         IPASS = 0
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
C     Restore KONTRL and check to see if the tests of error detection
C     passed.
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN,  9040)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN,  9050)
         ENDIF
      ENDIF
C
C     Print PASS/FAIL message.
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN, 9100)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN, 9110)
      RETURN
C
 9000 FORMAT ('1TEST OF SUBROUTINE LSEI')
 9010 FORMAT (/' LSEI PASSED TEST')
 9020 FORMAT (/' LSEI FAILED TEST'/' RELERR = ',1P,E20.6/' RELNRM = ',
     *        E20.6)
 9030 FORMAT (/ ' 2 ERROR MESSAGES EXPECTED')
 9040 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9050 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
 9100 FORMAT (/' ****************LSEI PASSED ALL TESTS***************')
 9110 FORMAT (/' ****************LSEI FAILED SOME TESTS**************')
      END
