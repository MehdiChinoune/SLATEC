*DECK DFZTST
      SUBROUTINE DFZTST (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  DFZTST
C***PURPOSE  Quick check for DFZERO.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (FZTEST-S, DFZTST-D)
C***AUTHOR  Boland, W. Robert, (LANL)
C***ROUTINES CALLED  D1MACH, DFZERO, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   920212  DATE WRITTEN
C***END PROLOGUE  DFZTST
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      INTEGER IFLAG, KONTRL
      DOUBLE PRECISION AE, B, C, PI, R, RE, TOL
      LOGICAL FATAL
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
C     .. External Subroutines ..
      EXTERNAL DFZERO, XERCLR, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, ATAN, DSIN, DTAN, MAX, SQRT
C***FIRST EXECUTABLE STATEMENT  DFZTST
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
      IPASS = 1
      PI = 4.0D0 *ATAN(1.0D0)
      RE = 1.0D-10
      AE = 1.0D-10
      TOL = MAX(1.0D-9,SQRT(D1MACH(4)))
C
C     Set up and solve example problem
C
      B = 0.1D0
      C = 4.0D0
      R = C - B
      CALL DFZERO (DSIN, B, C, R, RE, AE, IFLAG)
C
C     See if test was passed.
C
      IF (ABS(B-PI).LE.TOL .AND. ABS(C-PI).LE.TOL) THEN
        IF (KPRINT .GE. 3) WRITE (LUN, 9010) 'PASSED', B, C, IFLAG
      ELSE
        IPASS = 0
        IF (KPRINT .GE. 2) WRITE (LUN, 9010) 'FAILED', B, C, IFLAG
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
C
      IF (KPRINT .GE. 3) WRITE (LUN,9020)
      B = 1.0D0
C
C     IFLAG=3 (Singular point)
C
      C = 2.0D0
      R = 0.5D0*(B+C)
      CALL DFZERO (DTAN, B, C, B, RE, AE, IFLAG)
      IF (IFLAG .NE. 3) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9030) IFLAG, 2
      ENDIF
C
C     IFLAG=4 (No sign change)
C
      B = -3.0D0
      C = -0.1D0
      R = 0.5D0*(B+C)
      CALL DFZERO (DSIN, B, C, R, RE, AE, IFLAG)
      IF (IFLAG .NE. 4) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9030) IFLAG, 4
      ENDIF
C
      CALL XERCLR
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
        IF (KPRINT .GE. 2) THEN
          WRITE (LUN, 9040)
        ENDIF
      ELSE
        IF (KPRINT .GE. 3) THEN
          WRITE (LUN, 9050)
        ENDIF
      ENDIF
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN,9060)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN,9070)
      RETURN
 9000 FORMAT ('1' / ' DFZERO QUICK CHECK')
 9010 FORMAT (' Accuracy test ', A /
     +        ' Example problem results:  (answer = PI),  B =', F20.14,
     +        ' C =', F20.14 / ' IFLAG =', I2)
 9020 FORMAT (/ ' IFLAG 3 and 4 tests')
 9030 FORMAT (/' IFLAG test FAILED.  IFLAG =', I2, ', but should ',
     +        'have been', I2)
 9040 FORMAT (/ ' At least IFLAG test failed')
 9050 FORMAT (/ ' All IFLAG tests passed')
 9060 FORMAT (/' ***************DFZERO PASSED ALL TESTS**************')
 9070 FORMAT (/' ***************DFZERO FAILED SOME TESTS*************')
      END
