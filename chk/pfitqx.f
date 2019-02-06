*DECK PFITQX
      SUBROUTINE PFITQX (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  PFITQX
C***PURPOSE  Quick check for POLFIT, PCOEF and PVALUE.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PFITQX-S, DPFITT-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  CMPARE, PASS, PCOEF, POLFIT, PVALUE, R1MACH,
C                    XERCLR, XGETF, XSETF
C***COMMON BLOCKS    CHECK
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890921  Realigned order of variables in the COMMON block.
C           (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900911  Test problem changed and cosmetic changes to code.  (WRB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4) and modified the
C           FORMATs.  (RWC)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920214  Code restructured to test for all values of KPRINT and to
C           provide more PASS/FAIL information.  (WRB)
C***END PROLOGUE  PFITQX
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Scalars in Common ..
      REAL EPS, RP, SVEPS, TOL
      INTEGER IERP, IERR, NORD, NORDP
C     .. Arrays in Common ..
      REAL R(11)
C     .. Local Scalars ..
      REAL YFIT
      INTEGER I, ICNT, M, MAXORD
C     .. Local Arrays ..
      REAL A(97), TC(5), W(11), X(11), Y(11), YP(5)
      INTEGER ITEST(9)
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     .. External Subroutines ..
      EXTERNAL CMPARE, PASS, PCOEF, POLFIT, PVALUE
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C     .. Common blocks ..
      COMMON /CHECK/ EPS, R, RP, SVEPS, TOL, NORDP, NORD, IERP, IERR
C***FIRST EXECUTABLE STATEMENT  PFITQX
      IF (KPRINT .GE. 2) WRITE (LUN,FMT=9000)
C
C     Initialize variables for testing passage or failure of tests
C
      DO 100 I = 1,9
        ITEST(I) = 0
  100 CONTINUE
      ICNT = 0
      TOL = SQRT(R1MACH(4))
      M = 11
      DO 110 I = 1,M
        X(I) = I - 6
        Y(I) = X(I)**4
  110 CONTINUE
C
C     Test POLFIT
C     Input EPS is negative - specified level
C
      W(1) = -1.0E0
      EPS = -0.01E0
      SVEPS = EPS
      MAXORD = 8
      NORDP = 4
      RP = 625.0E0
      IERP = 1
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 130
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 130
      WRITE (LUN,FMT=9010)
      WRITE (LUN,FMT=9020)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 120
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  120 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Input EPS is negative - computed level
C
  130 EPS = -1.0E0
      SVEPS = EPS
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 150
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 150
      WRITE (LUN,FMT=9050)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 140
      WRITE (LUN,FMT=9060) MAXORD
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  140 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Input EPS is zero
C
  150 W(1) = -1.0E0
      EPS = 0.0E0
      SVEPS = EPS
      NORDP = 5
      MAXORD = 5
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 170
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 170
      WRITE (LUN,FMT=9070)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 160
      WRITE (LUN,FMT=9060) MAXORD
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  160 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Input EPS is positive
C
  170 IERP = 1
      NORDP = 4
      EPS = 75.0E0*R1MACH(4)
      SVEPS = EPS
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 190
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 190
      WRITE (LUN,FMT=9080)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 180
      WRITE (LUN,FMT=9060) MAXORD
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  180 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Improper input
C
  190 IERP = 2
      M = -2
C
C     Check for suppression of printing.
C
      CALL XGETF (KONTRL)
      IF (KPRINT .LE. 2) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
      CALL XERCLR
C
      IF (KPRINT .GE. 3) WRITE (LUN,9090)
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      ICNT = ICNT + 1
      IF (IERR .EQ. 2) THEN
        ITEST(ICNT) = 1
        IF (KPRINT .GE. 3) THEN
          WRITE (LUN, 9100) 'PASSED', IERR
        ENDIF
      ELSE
        IF (KPRINT .GE. 2) THEN
          WRITE (LUN, 9100) 'FAILED', IERR
        ENDIF
      ENDIF
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 210
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 210
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 200
C
C     Send message indicating passage or failure of test
C
  200 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
      CALL XERCLR
      CALL XSETF (KONTRL)
C
C     MAXORD too small to meet RMS error
C
  210 M = 11
      W(1) = -1.0E0
      EPS = 5.0E0*R1MACH(4)
      SVEPS = EPS
      RP = 553.0E0
      MAXORD = 2
      IERP = 3
      NORDP = 2
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 230
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 230
      WRITE (LUN,FMT=9110)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 220
      WRITE (LUN,FMT=9060) MAXORD
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  220 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     MAXORD too small to meet statistical test
C
  230 NORDP = 4
      IERP = 4
      RP = 625.0E0
      EPS = -0.01E0
      SVEPS = EPS
      MAXORD = 5
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
C
C     See if test passed
C
      CALL CMPARE (ICNT, ITEST)
C
C     Check for suppression of printing.
C
      IF (KPRINT .EQ. 0) GO TO 250
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 250
      WRITE (LUN,FMT=9120)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 240
      WRITE (LUN,FMT=9060) MAXORD
      WRITE (LUN,FMT=9030) SVEPS,NORDP,RP,IERP
      WRITE (LUN,FMT=9040) EPS,NORD,R(11),IERR
C
C     Send message indicating passage or failure of test
C
  240 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Test PCOEF
C
  250 MAXORD = 6
      EPS = 0.0E0
      SVEPS = EPS
      Y(6) = 1.0E0
      DO 260 I = 1,M
        W(I) = 1.0E0/(Y(I)**2)
  260 CONTINUE
      Y(6) = 0.0E0
      CALL POLFIT (M, X, Y, W, MAXORD, NORD, EPS, R, IERR, A)
      CALL PCOEF (4, 5.0E0, TC, A)
C
C     See if test passed
C
      ICNT = ICNT + 1
      IF (ABS(R(11)-TC(1)) .LE. TOL) ITEST(ICNT) = 1
C
C     Check for suppression of printing
C
      IF (KPRINT .EQ. 0) GO TO 280
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 280
      WRITE (LUN,FMT=9130)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 270
      WRITE (LUN,FMT=9140) R(11),TC(1)
C
C     Send message indicating passage or failure of test
C
  270 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Test PVALUE
C     Normal call
C
  280 CALL PVALUE (6, 0, X(8), YFIT, YP, A)
C
C     See if test passed
C
      ICNT = ICNT + 1
      IF (ABS(R(8)-YFIT) .LE. TOL) ITEST(ICNT) = 1
C
C     Check for suppression of printing
C
      IF (KPRINT .EQ. 0) GO TO 300
      IF (KPRINT.EQ.1 .AND. ITEST(ICNT).EQ.1) GO TO 300
      WRITE (LUN,FMT=9150)
      WRITE (LUN,FMT=9160)
      IF (KPRINT.LE.2 .AND. ITEST(ICNT).EQ.1) GO TO 290
      WRITE (LUN,FMT=9170) X(8),R(8),YFIT
C
C     Send message indicating passage or failure of test
C
  290 CALL PASS (LUN, ICNT, ITEST(ICNT))
C
C     Check to see if all tests passed
C
  300 IPASS = 1
      DO 310 I = 1,9
        IPASS = IPASS*ITEST(I)
  310 CONTINUE
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.3) WRITE (LUN,FMT=9180)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.2) WRITE (LUN,FMT=9190)
      RETURN
C
 9000 FORMAT ('1' / ' Test POLFIT, PCOEF and PVALUE')
 9010 FORMAT (' Exercise POLFIT')
 9020 FORMAT (' Input EPS is negative - specified significance level')
 9030 FORMAT (' Input EPS =  ', E15.8, '   correct order =  ', I3,
     +        '   R(1) = ', E15.8, '   IERR = ', I1)
 9040 FORMAT (' Output EPS = ', E15.8, '   computed order = ', I3,
     +        '   R(1) = ', E15.8, '   IERR = ', I1)
 9050 FORMAT (/ ' Input EPS is negative - computed significance level')
 9060 FORMAT (' Maximum order = ', I2)
 9070 FORMAT (/ ' Input EPS is zero')
 9080 FORMAT (/ ' Input EPS is positive')
 9090 FORMAT (/ ' Invalid input')
 9100 FORMAT (' POLFIT incorrect argument test ', A /
     +        ' IERR should be 2.  It is ', I4)
 9110 FORMAT (/ ' Cannot meet RMS error requirement')
 9120 FORMAT (/ ' Cannot satisfy statistical test')
 9130 FORMAT (/ ' Exercise PCOEF')
 9140 FORMAT (/ ' For C=1.0, correct coefficient = ', E15.8,
     +        '   computed = ', E15.8)
 9150 FORMAT (/ ' Exercise PVALUE')
 9160 FORMAT (' Normal execution')
 9170 FORMAT (' For X = ', F5.2, '   correct P(X) = ', E15.8,
     +        '    P(X) from PVALUE = ', E15.8)
 9180 FORMAT (/' ***************POLFIT PASSED ALL TESTS***************')
 9190 FORMAT (/' ***************POLFIT FAILED SOME TESTS**************')
      END
