*DECK CQRTST
      SUBROUTINE CQRTST (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  CQRTST
C***PURPOSE  Quick check for CPQR79.
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (RQRTST-S, CQRTST-C)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  CPQR79, NUMXER, PASS, R1MACH, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
C   911010  Code reworked and simplified.  (RWC and WRB)
C***END PROLOGUE  CQRTST
      INTEGER ITEST(2), ITMP(7)
      REAL WORK(144)
      COMPLEX COEFF1(9), COEFF2(2), COEFF3(2), ROOT(8), CHK1(8), CHK2
      LOGICAL FATAL
C
      DATA COEFF1 / (1.0,0.0), (-7.0,-2.0), (8.0,6.0), (28.0, 8.0),
     *              (-49.0,-24.0), (7.0,2.0), (-8.0,-6.0),
     *              (-28.0,-8.0), (48.0,24.0)/
      DATA COEFF2 / (1.0,1.0), (1.0,3.0) /
      DATA COEFF3 / (0.0,0.0), (1.0,3.0) /
      DATA CHK1 / (4.0,2.0), (3.0,0.0), (-2.0,0.0), (2.0,0.0),
     *            (0.0,-1.0), (-1.0,0.0), (0.0,1.0), (1.0,0.0) /
      DATA CHK2 / (-2.0,-1.0) /
C***FIRST EXECUTABLE STATEMENT  CQRTST
      IF (KPRINT .GE. 2) WRITE (LUN, 90000)
      TOL = SQRT(R1MACH(4))
      IPASS = 1
C
C     First test.
C
      CALL CPQR79 (8, COEFF1, ROOT, IERR, WORK)
C
C     Check to see if test passed.
C
      DO 10 I=1,7
         ITMP(I) = 0
   10 CONTINUE
C
C     Check for roots in any order.
C
      DO 30 I=1,7
         DO 20 J=1,7
            IF (ABS(ROOT(I)-CHK1(J)) .LE. TOL) THEN
               ITMP(J) = 1
               GOTO 30
            ENDIF
   20    CONTINUE
   30 CONTINUE
C
C     Check that we found all 7 roots.
C
      ITEST(1) = 1
      DO 40 I=1,7
         ITEST(1) = ITEST(1)*ITMP(I)
   40 CONTINUE
C
C     Print test results.
C
      IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.ITEST(1).EQ.0)) THEN
         WRITE (LUN, 90010)
         WRITE (LUN, 90020) (J,COEFF1(J), J=1,9)
         WRITE (LUN, 90030)
         WRITE (LUN, 90040) (J,ROOT(J), J=1,7)
      ENDIF
      IF (KPRINT .GE. 2) THEN
         CALL PASS (LUN, 1, ITEST(1))
      ENDIF
C
C     Set up next problem.
C
      CALL CPQR79 (1, COEFF2, ROOT, IERR, WORK)
C
C     Check to see if test passed.
C
      ITEST(2) = 1
      IF (ABS(ROOT(1)-CHK2) .GT. TOL) ITEST(2) = 0
C
C     Print test results for second test.
C
      IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.ITEST(1).EQ.0)) THEN
         WRITE (LUN, 90050)
         WRITE (LUN, 90010)
         WRITE (LUN, 90020) (J,COEFF2(J), J=1,2)
         WRITE (LUN, 90030)
         WRITE (LUN, 90040) (J,ROOT(J), J=1,1)
      ENDIF
      IF (KPRINT .GE. 2) THEN
         CALL PASS (LUN, 2, ITEST(2))
      ENDIF
C
C     Trigger 2 error conditions
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
      FATAL = .FALSE.
      CALL XERCLR
      IF (KPRINT .GE. 3) WRITE (LUN, 90060)
C
C     CALL CPQR79 with 0 degree polynomial.
C
      CALL CPQR79 (0, COEFF2, ROOT, IERR, WORK)
      IF (NUMXER(NERR) .NE. 3) THEN
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
C     CALL CPQR79 with zero leading coefficient.
C
      CALL CPQR79 (2, COEFF3, ROOT, IERR, WORK)
      IF (NUMXER(NERR) .NE. 2) THEN
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IPASS = 0
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 90070)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 90080)
         ENDIF
      ENDIF
C
C     See if all tests passed.
C
      IPASS = IPASS*ITEST(1)*ITEST(2)
C
      IF (IPASS.EQ.1 .AND. KPRINT.GT.1) WRITE (LUN,90100)
      IF (IPASS.EQ.0 .AND. KPRINT.NE.0) WRITE (LUN,90110)
      RETURN
C
90000 FORMAT ('1', /,' CPQR79 QUICK CHECK')
90010 FORMAT (/, ' CHECK REAL AND IMAGINARY PARTS OF ROOT' /
     *          ' COEFFICIENTS')
90020 FORMAT (/ (I6, 3X, 1P, 2E22.14))
90030 FORMAT (// 25X, 'TABLE of ROOTS' //
     *        '   ROOT         REAL  PART', 12X, 'IMAG  PART' /
     *        '  NUMBER', 8X, 2(' of  ZERO ', 12X))
90040 FORMAT (I6, 3X, 1P, 2E22.14)
90050 FORMAT (/, ' TEST SUBSEQUENT RELATED CALL')
90060 FORMAT (// ' TRIGGER 2 ERROR CONDITIONS' //)
90070 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/' **************CPQR79 PASSED ALL TESTS**************')
90110 FORMAT (/' **************CPQR79 FAILED SOME TESTS*************')
      END
