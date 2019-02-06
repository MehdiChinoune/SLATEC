*DECK DEG8CK
      SUBROUTINE DEG8CK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DEG8CK
C***PURPOSE  Quick check for DEXINT and DGAUS8.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (EG8CK-S, DEG8CK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C   DEG8CK is a quick check routine for DEXINT and DGAUS8.  Exponential
C   integrals from DEXINT are checked against quadratures from DGAUS8.
C
C***ROUTINES CALLED  D1MACH, DEXINT, DFEIN, DGAUS8
C***COMMON BLOCKS    DFEINX
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890718  Added check when testing error conditions.  (WRB)
C   890718  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   910708  Code revised to test error returns for all values of
C           KPRINT.  (WRB)
C   920206  Corrected argument list in CALL to DEXINT.  (WRB)
C***END PROLOGUE  DEG8CK
      COMMON /DFEINX/ X, A, FKM
      INTEGER I, ICASE, IE, IERR, II, IK, IPASS, IX, IY, K, KE, KK,
     *        KODE, KX, LUN, M, N, NM, NZ
      DOUBLE PRECISION A, ANS, ATOL, BB, EN, ER, EX, FKM, SIG, SUM,
     *                 TOL, T1, T2, X, XX, Y
      DOUBLE PRECISION D1MACH, DFEIN
      DIMENSION EN(4), Y(4), XX(5)
      LOGICAL FATAL
      EXTERNAL DFEIN
C***FIRST EXECUTABLE STATEMENT  DEG8CK
      IF (KPRINT .GE. 2) WRITE (LUN, 90000)
      IPASS=1
      TOL = SQRT(MAX(D1MACH(4),1.0D-18))
      DO 150 KODE=1,2
        IK = KODE - 1
        FKM = IK
        DO 140 N=1,25,8
          DO 130 M=1,4
            NM = N + M - 1
            DO 120 IX=1,25,8
              X = IX- 0.20D0
              CALL DEXINT(X, N, KODE, M, TOL, EN, NZ, IERR)
              KX = X+0.5D0
              IF (KX.EQ.0) KX = 1
              ICASE = 1
              A = N
              IF (KX.LE.N) GO TO 10
              ICASE = 2
              A = NM
              IF (KX.GE.NM) GO TO 10
              ICASE = 3
              A = KX
   10         CONTINUE
              SIG = 3.0D0/X
              T2 = 1.0D0
              SUM = 0.0D0
   20         CONTINUE
              T1 = T2
              T2 = T2 + SIG
              ATOL = TOL
              CALL DGAUS8(DFEIN, T1, T2, ATOL, ANS, IERR)
              SUM = SUM + ANS
              IF (ABS(ANS).LT.ABS(SUM)*TOL) GO TO 30
              GO TO 20
   30         CONTINUE
              EX = 1.0D0
              IF (KODE.EQ.1) EX = EXP(-X)
              BB = A
              IF (ICASE.NE.3) GO TO 40
              IY = KX - N + 1
              Y(IY) = SUM
              KE = M - IY
              IE = IY - 1
              KK = IY
              II = IY
              GO TO 60
   40         CONTINUE
              IF (ICASE.NE.2) GO TO 50
              Y(M) = SUM
              IF (M.EQ.1) GO TO 100
              IE = M - 1
              II = M
              GO TO 80
   50         CONTINUE
              Y(1) = SUM
              IF (M.EQ.1) GO TO 100
              KE = M - 1
              KK = 1
   60         CONTINUE
C
C             Forward recur
C
              DO 70 K=1,KE
                Y(KK+1) = (EX-X*Y(KK))/BB
                BB = BB + 1.0D0
                KK = KK + 1
   70         CONTINUE
              IF (ICASE.NE.3) GO TO 100
   80         BB = A - 1.0D0
C
C             Backward recur
C
              DO 90 I=1,IE
                Y(II-1) = (EX-BB*Y(II))/X
                BB = BB - 1.0D0
                II = II - 1
   90         CONTINUE
  100         CONTINUE
              DO 110 I=1,M
                ER = ABS((Y(I)-EN(I))/Y(I))
                IF (ER .GT. TOL) THEN
                   WRITE (LUN,90010)
                   IPASS = 0
                   GO TO 160
                ENDIF
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
C
C     Trigger 6 error conditions.
C
  160 FATAL = .FALSE.
C
      IF (KPRINT .GE. 3) WRITE (LUN, 90020)
      XX(1) = 1.0D0
      XX(2) = 1.0D0
      XX(3) = 1.0D0
      XX(4) = 1.0D0
      XX(5) = 0.01D0
      DO 170 I=1,5
        XX(I) = -XX(I)
        K = XX(2)
        N = XX(3)
        M = XX(4)
        CALL DEXINT (XX(I), N, K, M, XX(5), EN, NZ, IERR)
        IF (IERR .NE. 1) THEN
           IPASS = 0
           FATAL = .TRUE.
           WRITE (LUN, 90030) I
        ENDIF
        XX(I) = -XX(I)
  170 CONTINUE
      X = 0.0D0
      TOL = 1.0D-2
      CALL DEXINT (X, 1, 1, 1, TOL, EN, NZ, IERR)
      IF (IERR .NE. 1) THEN
         IPASS = 0
         FATAL = .TRUE.
         WRITE (LUN, 90040)
      ENDIF
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 90070)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 90080)
         ENDIF
      ENDIF
C
      IF(IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN, 90100)
      IF(IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN, 90110)
      RETURN
C
90000 FORMAT ('1' / ' QUICK CHECK FOR DEXINT AND DGAUS8' /)
90010 FORMAT (// ' ERROR IN DEG8CK COMPARISON TEST' /)
90020 FORMAT (/ ' TRIGGER 6 ERROR CONDITIONS' /)
90030 FORMAT (' Error occurred with DO index I =', I2)
90040 FORMAT (' Error occurred with X = 0.0')
90070 FORMAT (/' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/ ' *********DEXINT AND DGAUS8 PASSED ALL TESTS*********')
90110 FORMAT (/ ' *********DEXINT OR DGAUS8 FAILED SOME TESTS*********')
      END
