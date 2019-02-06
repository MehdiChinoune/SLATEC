*DECK ITEST
      SUBROUTINE ITEST (LEN, ICOMP, ITRUE, KPRINT)
C***BEGIN PROLOGUE  ITEST
C***PURPOSE  Compare arrays ICOMP and ITRUE.
C***LIBRARY   SLATEC
C***TYPE      INTEGER (ITEST-I)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Lawson, C. L., (JPL)
C***DESCRIPTION
C
C   This subroutine compares the arrays ICOMP and ITRUE of length LEN
C   for equality.  In the case of an unequal compare, appropriate
C   messages are written.
C
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    COMBLA
C***REVISION HISTORY  (YYMMDD)
C   741210  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920211  Code restructured and information added to the DESCRIPTION
C           section.  (WRB)
C***END PROLOGUE  ITEST
      INTEGER ICOMP(*), ITRUE(*)
      LOGICAL PASS
      COMMON /COMBLA/ NPRINT, ICASE, N, INCX, INCY, MODE, PASS
C***FIRST EXECUTABLE STATEMENT  ITEST
      DO 100 I = 1,LEN
        IF (ICOMP(I) .NE. ITRUE(I)) THEN
C
C         Here ICOMP(I) is not equal to ITRUE(I).
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
          IF (KPRINT .GE. 3) THEN
            ID = ICOMP(I) - ITRUE(I)
            WRITE (NPRINT,9020) ICASE, N, INCX, INCY, MODE, I, ICOMP(I),
     +                          ITRUE(I), ID
          ENDIF
        ENDIF
  100 CONTINUE
      RETURN
 9000 FORMAT ('+', 39X, 'FAIL')
 9010 FORMAT ('0CASE  N INCX INCY MODE  I', 29X, 'COMP(I)', 29X,
     +        'TRUE(I)', 2X, 'DIFFERENCE' / 1X)
 9020 FORMAT (1X, I4, I3, 3I5, I3, 2I36, I12)
      END
