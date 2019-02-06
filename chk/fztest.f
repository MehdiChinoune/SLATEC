*DECK FZTEST
      SUBROUTINE FZTEST (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  FZTEST
C***PURPOSE  Quick check for FZERO.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FZTEST-S, DFZTST-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  FZERO, R1MACH, XERCLR, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
C   910501  Added PURPOSE and TYPE records.  (WRB)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920212  Code completely restructured to test IFLAG for all values
C           of KPRINT.  (WRB)
C***END PROLOGUE  FZTEST
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      INTEGER IFLAG, KONTRL
      REAL AE, B, C, PI, R, RE, TOL
      LOGICAL FATAL
C     .. External Functions ..
      REAL R1MACH
      EXTERNAL R1MACH
C     .. External Subroutines ..
      EXTERNAL FZERO, XERCLR, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, ATAN, MAX, SIN, SQRT, TAN
C***FIRST EXECUTABLE STATEMENT  FZTEST
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
      IPASS = 1
      PI = 4.0E0 *ATAN(1.0E0)
      RE = 1.0E-6
      AE = 1.0E-6
      TOL = MAX(1.0E-5,SQRT(R1MACH(4)))
C
C     Set up and solve example problem
C
      B = 0.1E0
      C = 4.0E0
      R = C - B
      CALL FZERO (SIN, B, C, R, RE, AE, IFLAG)
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
      B = 1.0E0
C
C     IFLAG=3 (Singular point)
C
      C = 2.0E0
      R = 0.5E0*(B+C)
      CALL FZERO (TAN, B, C, B, RE, AE, IFLAG)
      IF (IFLAG .NE. 3) THEN
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9030) IFLAG, 2
      ENDIF
C
C     IFLAG=4 (No sign change)
C
      B = -3.0E0
      C = -0.1E0
      R = 0.5E0*(B+C)
      CALL FZERO (SIN, B, C, R, RE, AE, IFLAG)
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
 9000 FORMAT ('1' / ' FZERO QUICK CHECK')
 9010 FORMAT (' Accuracy test ', A /
     +        ' Example problem results:  (answer = PI),  B =', F20.14,
     +        ' C =', F20.14 / ' IFLAG =', I2)
 9020 FORMAT (/ ' IFLAG 3 and 4 tests')
 9030 FORMAT (/' IFLAG test FAILED.  IFLAG =', I2, ', but should ',
     +        'have been', I2)
 9040 FORMAT (/ ' At least IFLAG test failed')
 9050 FORMAT (/ ' All IFLAG tests passed')
 9060 FORMAT (/' ***************FZERO PASSED ALL TESTS**************')
 9070 FORMAT (/' ***************FZERO FAILED SOME TESTS*************')
      END
