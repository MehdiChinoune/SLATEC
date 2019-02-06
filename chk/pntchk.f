*DECK PNTCHK
      SUBROUTINE PNTCHK (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  PNTCHK
C***PURPOSE  Quick check for POLINT, POLCOF and POLYVL
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PNTCHK-S, DPNTCK-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  NUMXER, POLCOF, POLINT, POLYVL, R1MACH, XERCLR,
C                    XGETF, XSETF
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
C   910708  Minor modifications in use of KPRINT.  (WRB)
C   920212  Code completely restructured to test errors for all values
C           of KPRINT.  (WRB)
C***END PROLOGUE  PNTCHK
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL TOL, YF
      INTEGER I, IERR, KONTRL, N, NERR
      LOGICAL FATAL
C     .. Local Arrays ..
      REAL C(6), D(6), DCHK(6), W(12), X(6), XCHK(6), Y(6)
C     .. External Functions ..
      REAL R1MACH
      INTEGER NUMXER
      EXTERNAL R1MACH, NUMXER
C     .. External Subroutines ..
      EXTERNAL POLCOF, POLINT, POLYVL, XERCLR, XGETF, XSETF
C     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
C     .. Data statements ..
      DATA X / 1.0E0, 2.0E0, 3.0E0, -1.0E0, -2.0E0, -3.0E0 /
      DATA Y / 0.0E0, 9.0E0, 64.0E0, 0.0E0, 9.0E0, 64.0E0 /
      DATA XCHK / 1.0E0, 0.0E0, -2.0E0, 0.0E0, 1.0E0, 0.0E0 /
      DATA DCHK / 1.0E0, 0.0E0, -4.0E0, 0.0E0, 24.0E0, 0.0E0 /
C***FIRST EXECUTABLE STATEMENT  PNTCHK
      IF (KPRINT .GE. 2) WRITE (LUN,9000)
C
C     Initialize variables for tests.
C
      TOL = SQRT(R1MACH(4))
      IPASS = 1
      N = 6
C
C     Set up polynomial test.
C
      CALL POLINT (N, X, Y, C)
      CALL POLCOF (0.0E0, N, X, C, D, W)
C
C     Check to see if POLCOF test passed.
C
      FATAL = .FALSE.
      DO 110 I = 1,N
        IF (ABS(D(I)-XCHK(I)) .GT. TOL) THEN
          IPASS = 0
          FATAL = .TRUE.
        ENDIF
  110 CONTINUE
      IF (FATAL) THEN
        IF (KPRINT .GE. 2) WRITE (LUN, 9010) 'FAILED', (D(I), I = 1,N)
      ELSE
        IF (KPRINT .GE. 3) WRITE (LUN, 9010) 'PASSED', (D(I), I = 1,N)
      ENDIF
C
C     Test POLYVL.
C
      CALL POLYVL (5, 0.0E0, YF, D, N, X, C, W, IERR)
      IF (ABS(DCHK(1)-YF) .LE. TOL) THEN
        IF (KPRINT .GE. 3) WRITE (LUN, 9020) 'PASSED', YF,(D(I),I=1,5)
      ELSE
        IPASS = 0
        IF (KPRINT .GE. 2) WRITE (LUN, 9020) 'FAILED', YF,(D(I),I=1,5)
      ENDIF
C
      FATAL = .FALSE.
      DO 120 I = 1,5
        IF (ABS(DCHK(I+1)-D(I)) .GT. TOL) THEN
          IPASS = 0
          FATAL = .TRUE.
        ENDIF
  120 CONTINUE
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
      IF (KPRINT .GE. 3) WRITE (LUN,9030)
      CALL POLINT (0, X, Y, C)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
      CALL XERCLR
C
      X(1) = -1.0E0
      CALL POLINT (N, X, Y, C)
      IF (NUMXER(NERR) .NE. 2) THEN
        IPASS = 0
        FATAL = .TRUE.
      ENDIF
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
      IF (IPASS.EQ.1 .AND. KPRINT.GE.2) WRITE (LUN,9080)
      IF (IPASS.EQ.0 .AND. KPRINT.GE.1) WRITE (LUN,9090)
      RETURN
C
 9000 FORMAT ('1' / ' Test POLINT, POLCOF and POLYVL')
 9010 FORMAT (/ 'POLCOF ', A, ' test' /
     +        ' Taylor coefficients for the quintic should be' /
     +        6X, '1.000', 5X, '0.000', 4X, '-2.000', 5X, '0.000', 5X,
     +        '1.000', 5X, '0.000' /
     +        ' Taylor coefficients from POLCOF are' / 1X, 6F10.3 /)
 9020 FORMAT (' Derivative test ', A /
     +        ' The derivatives of the polynomial at zero as ',
     +        'computed by POLYVL are' / 1X, 6F10.3 /)
 9030 FORMAT (/' 2 Error messages expected')
 9040 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
 9050 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
 9080 FORMAT (/' ****************POLINT PASSED ALL TESTS**************')
 9090 FORMAT (/' ***************POLINT FAILED SOME TESTS**************')
      END
