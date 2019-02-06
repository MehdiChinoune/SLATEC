*DECK STEST
      SUBROUTINE STEST (LEN, SCOMP, STRUE, SSIZE, SFAC, KPRINT)
C***BEGIN PROLOGUE  STEST
C***PURPOSE  Compare arrays SCOMP and STRUE.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (STEST-S, DTEST-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Lawson, C. L., (JPL)
C***DESCRIPTION
C
C   This subroutine compares arrays SCOMP and STRUE of length LEN to
C   see if the term by term differences, multiplied by SFAC, are
C   negligible.  In the case of a significant difference, appropriate
C   messages are written.
C
C***ROUTINES CALLED  R1MACH
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   741210  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900820  Modified IF test to use function DIFF and made cosmetic
C           changes to routine.  (WRB)
C   901005  Removed usage of DIFF in favour of R1MACH.  (RWC)
C   910501  Added TYPE record.  (WRB)
C   920211  Code restructured and information added to the DESCRIPTION
C           section.  (WRB)
C***END PROLOGUE  STEST
      REAL SCOMP(*), STRUE(*), SSIZE(*), SFAC, SD, RELEPS, R1MACH
      LOGICAL PASS
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
      SAVE RELEPS
      DATA RELEPS /0.0E0/
C***FIRST EXECUTABLE STATEMENT  STEST
      IF (RELEPS .EQ. 0.0E0) RELEPS = R1MACH(4)
      DO 100 I = 1,LEN
        SD = ABS(SCOMP(I)-STRUE(I))
        IF (SFAC*SD .GT. ABS(SSIZE(I))*RELEPS) THEN
C
C         Here SCOMP(I) is not close to STRUE(I).
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
     +                       MODE, I, SCOMP(I), STRUE(I), SD, SSIZE(I)
        ENDIF
  100 CONTINUE
      RETURN
 9000 FORMAT ('+', 39X, 'FAIL')
 9010 FORMAT ('0CASE  N INCX INCY MODE  I', 29X, 'COMP(I)', 29X,
     +        'TRUE(I)', 2X, 'DIFFERENCE', 5X, 'SIZE(I)' / 1X)
 9020 FORMAT (1X, I4, I3, 3I5, I3, 2E36.8, 2E12.4)
      END
