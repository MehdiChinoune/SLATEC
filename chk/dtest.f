*DECK DTEST
      SUBROUTINE DTEST (LEN, DCOMP, DTRUE, DSIZE, DFAC, KPRINT)
C***BEGIN PROLOGUE  DTEST
C***PURPOSE  Compare arrays DCOMP and DTRUE.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (STEST-S, DTEST-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Lawson, C. L., (JPL)
C***DESCRIPTION
C
C   This subroutine compares arrays DCOMP and DTRUE of length LEN to
C   see if the term by term differences, multiplied by DFAC, are
C   negligible.  In the case of a significant difference, appropriate
C   messages are written.
C
C***ROUTINES CALLED  D1MACH
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   741210  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900820  Modified IF test to use function DDIFF and made cosmetic
C           changes to routine.  (WRB)
C   901005  Removed usage of DDIFF in favour of D1MACH.  (RWC)
C   910501  Added TYPE record.  (WRB)
C   920211  Code restructured and information added to the DESCRIPTION
C           section.  (WRB)
C***END PROLOGUE  DTEST
      DOUBLE PRECISION DCOMP(*), DTRUE(*), DSIZE(*), DFAC, DD,
     +       RELEPS, D1MACH
      LOGICAL PASS
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
      SAVE RELEPS
      DATA RELEPS /0.0D0/
C***FIRST EXECUTABLE STATEMENT  DTEST
      IF (RELEPS .EQ. 0.0D0) RELEPS = D1MACH(4)
      DO 100 I = 1,LEN
        DD = ABS(DCOMP(I)-DTRUE(I))
        IF (DFAC*DD .GT. ABS(DSIZE(I))*RELEPS) THEN
C
C         Here DCOMP(I) is not close to DTRUE(I).
C
          IF (PASS) THEN
C
C           Print FAIL message and header.
C
            PASS = .FALSE.
            IF (KPRINT .GE. 3) THEN
              WRITE (NPRINT,9000)
              WRITE (NPRINT,9010)
            ENDIF
          ENDIF
          IF (KPRINT .GE. 3) WRITE (NPRINT,9020) ICASE, N, INCX, INCY,
     +                       MODE, I, DCOMP(I), DTRUE(I), DD, DSIZE(I)
        ENDIF
  100 CONTINUE
      RETURN
 9000 FORMAT ('+', 39X, 'FAIL')
 9010 FORMAT ('0CASE  N INCX INCY MODE  I', 29X, 'COMP(I)', 29X,
     +        'TRUE(I)', 2X, 'DIFFERENCE', 5X, 'SIZE(I)' / 1X)
 9020 FORMAT (1X, I4, I3, 3I5, I3, 2D36.18, 2D12.4)
      END
