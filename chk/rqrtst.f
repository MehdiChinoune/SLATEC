*DECK RQRTST
      SUBROUTINE RQRTST (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  RQRTST
C***PURPOSE  Quick check for RPQR79.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (RQRTST-S, CQRTST-C)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  NUMXER, PASS, R1MACH, RPQR79, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901010  Restructured using IF-THEN-ELSE-ENDIF, cleaned up FORMATs
C           and changed TOL from sqrt R1MACH(3) to sqrt R1MACH(4) for
C           the IBM 370 mainframes.  (RWC)
C   911010  Code reworked and simplified.  (RWC and WRB)
C***END PROLOGUE  RQRTST
      INTEGER ITMP(7)
      COMPLEX ROOT(7), CHK(7)
      DIMENSION WORK(63)
      REAL COEF(8)
      LOGICAL FATAL
C
      DATA CHK / ( 1.4142135623731,  1.4142135623731),
     *           ( 1.4142135623731, -1.4142135623731),
     *           (0.0, 2.0), (0.0, -2.0), (-2.0, 0.0),
     *           (-1.4142135623731,  1.4142135623731),
     *           (-1.4142135623731, -1.4142135623731) /
C***FIRST EXECUTABLE STATEMENT  RQRTST
      IF (KPRINT .GE. 2) WRITE (LUN, 90000)
      TOL = SQRT(R1MACH(4))
      IPASS = 1
C
C     Initialize variables for testing.
C
      BETA = 0.0078125
      DO 20 J=1,8
         COEF(J) = BETA
         BETA = 2.0*BETA
   20 CONTINUE
C
      CALL RPQR79 (7, COEF, ROOT, IERR, WORK)
C
C     Check to see if test passed.
C
      DO 10 I=1,7
         ITMP(I) = 0
   10 CONTINUE
C
C     Check for roots in any order.
C
      DO 40 I=1,7
         DO 30 J=1,7
            IF (ABS(ROOT(I)-CHK(J)) .LE. TOL) THEN
               ITMP(J) = 1
               GO TO 40
            ENDIF
   30    CONTINUE
   40 CONTINUE
C
C     Check that we found all 7 roots.
C
      IPASS = 1
      DO 50 I=1,7
         IPASS = IPASS*ITMP(I)
   50 CONTINUE
C
C     Print test results.
C
      IF (KPRINT.GE.3 .OR. (KPRINT.GE.2.AND.IPASS.EQ.0)) THEN
         WRITE (LUN, 90010)
         WRITE (LUN, 90020) (J,COEF(J), J=1,8)
         WRITE (LUN, 90030)
         WRITE (LUN, 90040) (J,ROOT(J), J=1,7)
      ENDIF
      IF (KPRINT .GE. 2) THEN
         CALL PASS (LUN, 1, IPASS)
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
C     CALL RPQR79 with 0 degree polynomial.
C
      CALL RPQR79 (0, COEF, ROOT, IERR, WORK)
      IF (NUMXER(NERR) .NE. 3) THEN
         FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
C     CALL RPQR79 with zero leading coefficient.
C
      COEF(1) = 0.0
      CALL RPQR79 (2, COEF, ROOT, IERR, WORK)
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
      IF (IPASS.EQ.1 .AND. KPRINT.GT.1) WRITE (LUN,90100)
      IF (IPASS.EQ.0 .AND. KPRINT.NE.0) WRITE (LUN,90110)
      RETURN
C
90000 FORMAT ('1', /,' RPQR79 QUICK CHECK')
90010 FORMAT (/, ' CHECK REAL AND IMAGINARY PARTS OF ROOT' /
     *          ' COEFFICIENTS')
90020 FORMAT (/ (I6, 3X, 1P, E22.14))
90030 FORMAT (// 25X, 'TABLE of ROOTS' //
     *        '   ROOT         REAL  PART', 12X, 'IMAG  PART' /
     *        '  NUMBER', 8X, 2(' of  ZERO ', 12X))
90040 FORMAT (I6, 3X, 1P, 2E22.14)
90060 FORMAT (// ' TRIGGER 2 ERROR CONDITIONS' //)
90070 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/' **************RPQR79 PASSED ALL TESTS**************')
90110 FORMAT (/' **************RPQR79 FAILED SOME TESTS*************')
      END
