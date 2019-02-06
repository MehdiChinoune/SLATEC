*DECK QG8TST
      SUBROUTINE QG8TST (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QG8TST
C***PURPOSE  Quick check for GAUS8.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QG8TST-S, DQG8TS-D)
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  FQD1, FQD2, GAUS8, R1MACH, XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
C   910501  Added PURPOSE and TYPE records.  (WRB)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920213  Code restructured to test GAUS8 for all values of KPRINT,
C           second accuracy test added and testing of error returns
C           revised.  (WRB)
C***END PROLOGUE  QG8TST
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      INTEGER IERR, KONTRL
      REAL A, ANS, B, COR, ERR, REQ, TOL
      LOGICAL FATAL
C     .. External Functions ..
      REAL FQD1, FQD2, R1MACH
      EXTERNAL FQD1, FQD2, R1MACH
C     .. External Subroutines ..
      EXTERNAL GAUS8, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, ATAN, EXP, SQRT
C***FIRST EXECUTABLE STATEMENT  QG8TST
      IF (KPRINT .GE. 2) WRITE (LUN,FMT=9000)
C
C     Initialize variables for testing.
C
      TOL = SQRT(R1MACH(4))
      IPASS = 1
C
C     First accuracy test.
C
      A = 1.0E0
      B = 4.0E0
      ERR = TOL/100.0E0
      CALL GAUS8 (FQD1, A, B, ERR, ANS, IERR)
      COR = 2.0E0
      IF (ABS(ANS-COR).LE.TOL .AND. IERR.EQ.1) THEN
        IF (KPRINT .GE. 3)
     +    WRITE (LUN, 9010) 'PASSED', A, B, ANS, COR, ERR, IERR
      ELSE
        IPASS = 0
        IF (KPRINT .GE. 2)
     +    WRITE (LUN, 9010) 'FAILED', A, B, ANS, COR, ERR, IERR
      ENDIF
C
C     Second accuracy test.
C
      A = 0.0E0
      B = 4.0E0*ATAN(1.0E0)
      ERR = TOL/100.0E0
      CALL GAUS8 (FQD2, A, B, ERR, ANS, IERR)
      COR = (EXP(B)-1.0E0)/101.0E0
      IF (ABS(ANS-COR).LE.TOL .AND. IERR.EQ.1) THEN
        IF (KPRINT .GE. 3)
     +    WRITE (LUN, 9010) 'PASSED', A, B, ANS, COR, ERR, IERR
      ELSE
        IPASS = 0
        IF (KPRINT .GE. 2)
     +    WRITE (LUN, 9010) 'FAILED', A, B, ANS, COR, ERR, IERR
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
C
      IF (KPRINT .GE. 3) WRITE (LUN,FMT=9030)
C
C     Test with a discontinuous integrand and a tight error tolerance.
C
      A = 0.0E0
      B = 1.0E0
      COR = 2.0E0
      ERR = 100.0E0*R1MACH(4)
      REQ = ERR
      CALL GAUS8 (FQD1, A, B, ERR, ANS, IERR)
C
C     See if test passed.
C
      IF (IERR .EQ. 2) THEN
        IF (KPRINT .GE. 3)
     +    WRITE (LUN,FMT=9040) 'PASSED', REQ, ANS, IERR, ERR, COR
      ELSE
        IF (KPRINT .GE. 2)
     +    WRITE (LUN,FMT=9040) 'PASSED', REQ, ANS, IERR, ERR, COR
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
C
C     Test GAUS8 with A and B nearly equal.
C
      A = 2.0E0
      B = A*(1.0E0+R1MACH(4))
      COR = 0.0E0
      ERR = TOL
C
      CALL GAUS8 (FQD1, A, B, ERR, ANS, IERR)
C
C     Check to see if test passed.
C
      IF (IERR.EQ.-1 .AND. ANS.EQ.0.0E0) THEN
        IF (KPRINT .GE. 3) WRITE (LUN,9050) 'PASSED'
      ELSE
        IPASS = 0
        FATAL = .TRUE.
        IF (KPRINT .GE. 2) WRITE (LUN,9050) 'FAILED'
      ENDIF
C
      CALL XSETF (KONTRL)
      IF (FATAL) THEN
         IF (KPRINT .GE. 2) THEN
            WRITE (LUN, 9060)
         ENDIF
      ELSE
         IF (KPRINT .GE. 3) THEN
            WRITE (LUN, 9070)
         ENDIF
      ENDIF
C
      IF (IPASS.EQ.1 .AND. KPRINT.GE.3) WRITE (LUN,FMT=9080)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.2) WRITE (LUN,FMT=9090)
      RETURN
C
 9000 FORMAT ('1' / ' GAUS8 Quick Check')
 9010 FORMAT (/ ' Accuracy test of GAUS8 ', A /
     +        ' A = ', F10.5, '   B = ', F10.5 /
     +        ' Computed result = ', E14.7, '   Exact result = ',
     +        E14.7 /
     +        ' Tolerance = ', E14.7, '   IERR = ', I2 /)
 9030 FORMAT (/ ' Test error returns' /
     +        ' 2 error messages expected' /)
 9040 FORMAT (' Test of GAUS8 ', A /
     +        ' REQ =', E10.2, 5X, 'ANS =', E20.13, 5X, 'IERR =', I2,
     +        5X, 'should be 2' /
     +        ' ERR =', E10.2, ' CORRECT =' ,E20.13 /)
 9050 FORMAT (' Test of A and B nearly equal ', A)
 9060 FORMAT (/ ' At least one incorrect argument test FAILED')
 9070 FORMAT (/ ' All incorrect argument tests PASSED')
 9080 FORMAT (/,' ***************GAUS8 PASSED ALL TESTS***************')
 9090 FORMAT (/,' ***************GAUS8 FAILED SOME TESTS**************')
      END
