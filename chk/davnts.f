*DECK DAVNTS
      SUBROUTINE DAVNTS (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DAVNTS
C***PURPOSE  Quick check for DAVINT.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (AVNTST-S, DAVNTS-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  D1MACH, DAVINT, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of D1MACH(3) to D1MACH(4).  (RWC)
C   910501  Added PURPOSE and TYPE records.  (WRB)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920210  Code restructured and revised to test error returns for all
C           values of KPRINT.  (WRB)
C***END PROLOGUE  DAVNTS
      DOUBLE PRECISION D1MACH
      INTEGER I, IERR, IPASS, KPRINT, LUN, N
      DOUBLE PRECISION A, ANS, B, DEL, RN1, SQB, TOL, TOL1, X(501),
     +                 XINT, Y(501)
      LOGICAL FATAL
C***FIRST EXECUTABLE STATEMENT  DAVNTS
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
      IPASS = 1
      TOL = MAX(.0001D0,SQRT(D1MACH(4)))
      TOL1 = 1.0D-2*TOL
C
C     Perform first accuracy test.
C
      A = 0.0D0
      B = 5.0D0
      XINT = EXP(5.0D0) - 1.0D0
      N = 500
      RN1 = N - 1
      SQB = SQRT(B)
      DEL = 0.4D0*(B-A)/(N-1)
      DO 100 I = 1,N
        X(I) = SQB*SQRT(A+(I-1)*(B-A)/RN1) + DEL
        Y(I) = EXP(X(I))
  100 CONTINUE
      CALL DAVINT (X, Y, N, A, B, ANS, IERR)
C
C     See if test was passed.
C
      IF (ABS(ANS-XINT) .GT. TOL) THEN
        IPASS = 0
        IF (KPRINT .GE. 3) WRITE (LUN,9010) IERR, ANS, XINT
      ENDIF
C
C     Perform second accuracy test.
C
      X(1) = 0.0D0
      X(2) = 5.0D0
      Y(1) = 1.0D0
      Y(2) = 0.5D0
      A = -0.5D0
      B = 0.5D0
      XINT = 1.0D0
      CALL DAVINT (X, Y, 2, A, B, ANS, IERR)
C
C     See if test was passed.
C
      IF (ABS(ANS-XINT) .GT. TOL1) THEN
        IPASS = 0
        IF (KPRINT .GE. 3) WRITE (LUN,9010) IERR, ANS, XINT
      ENDIF
C
C     Send message indicating passage or failure of tests.
C
      IF (KPRINT .GE. 2) THEN
        IF (IPASS .EQ. 1) THEN
           IF (KPRINT .GE. 3) WRITE (LUN,9020)
        ELSE
           WRITE (LUN,9030)
        ENDIF
      ENDIF
C
C     Test error returns.
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
      IF (KPRINT .GE. 3) THEN
        WRITE (LUN,9040)
      ENDIF
      DO 110 I = 1,20
        X(I) = (I-1)/19.0D0 - 0.01D0
        IF (I .NE. 1) Y(I) = X(I)/(EXP(X(I))-1.0)
  110 CONTINUE
C
C     Test IERR = 1 error return.
C
      Y(1) = 1.0D0
      CALL DAVINT (X, Y, 20, 0.0D0, 1.0D0, ANS, IERR)
      IF (IERR .NE. 1) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9060) IERR, 1
      ENDIF
      CALL XERCLR
C
C     Test IERR = 2 error return.
C
      CALL DAVINT (X, Y, 20, 1.0D0, 0.0D0, ANS, IERR)
      IF (IERR .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9060) IERR, 2
      ENDIF
      IF (ANS .NE. 0.0D0) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9070)
      ENDIF
      CALL XERCLR
C
C     Test IERR = 5 error return.
C
      CALL DAVINT (X, Y, 1, 0.0D0, 1.0D0, ANS, IERR)
      IF (IERR .NE. 5) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9060) IERR, 5
      ENDIF
      IF (ANS .NE. 0.0D0) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9070)
      ENDIF
      CALL XERCLR
C
C     Test IERR = 4 error return.
C
      X(1) = 1.0D0/19.0D0
      X(2) = 0.0D0
      CALL DAVINT (X, Y, 20, 0.0D0, 1.0D0, ANS, IERR)
      IF (IERR .NE. 4) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9060) IERR, 4
      ENDIF
      IF (ANS .NE. 0.0D0) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9070)
      ENDIF
      CALL XERCLR
C
C     Test IERR = 3 error return.
C
      X(1) = 0.0D0
      X(2) = 1.0D0/19.0D0
      CALL DAVINT (X, Y, 20, 0.0D0, .01D0, ANS, IERR)
      IF (IERR .NE. 3) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9060) IERR, 3
      ENDIF
      IF (ANS .NE. 0.0D0) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 3) WRITE (LUN,9070)
      ENDIF
      CALL XERCLR
C
C     Reset XERMSG control variables and write summary.
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 9080)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 9090)
         ENDIF
      ENDIF
C
C     Write PASS/FAIL message.
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.3) WRITE (LUN,9100)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.2) WRITE (LUN,9110)
      RETURN
 9000 FORMAT ('1' / ' DAVINT Quick Check')
 9010 FORMAT (/' FAILED ACCURACY TEST' /
     +        ' IERR=', I2, 5X, 'COMPUTED ANS=', E20.11 / 14X,
     +        'CORRECT ANS=', D20.11, 5X, 'REQUESTED ERR=', D10.2)
 9020 FORMAT (/ ' DAVINT passed both accuracy tests.')
 9030 FORMAT (/ ' DAVINT failed at least one accuracy test.')
 9040 FORMAT (/ ' Test error returns from DAVINT' /
     +        ' 4 error messages expected' /)
 9060 FORMAT (/' IERR =', I2, ' and it should =', I2 /)
 9070 FORMAT (1X, 'ANS .NE. 0')
 9080 FORMAT (/ ' At least one incorrect argument test FAILED')
 9090 FORMAT (/ ' All incorrect argument tests PASSED')
 9100 FORMAT (/' ***************DAVINT PASSED ALL TESTS***************')
 9110 FORMAT (/' ***************DAVINT FAILED SOME TESTS**************')
      END
