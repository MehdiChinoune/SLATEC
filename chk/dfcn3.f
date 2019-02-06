*DECK DFCN3
      SUBROUTINE DFCN3 (IFLAG, M, N, X, FVEC, FJROW, NROW)
C***BEGIN PROLOGUE  DFCN3
C***PURPOSE  Subsidiary to DNLS1Q.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (FCN3-S, DFCN3-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Subroutine to evaluate the Jacobian, one row at a time, for
C   test problem used in quick check of DNLS1E.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  TYPE and declarations sections added and code polished.
C           (WRB)
C***END PROLOGUE  DFCN3
C     .. Scalar Arguments ..
      INTEGER IFLAG, M, N, NROW
C     .. Array Arguments ..
      DOUBLE PRECISION FJROW(*), FVEC(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP, TWO
      INTEGER I
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     .. Data statements ..
      DATA TWO /2.0D0/
C***FIRST EXECUTABLE STATEMENT  DFCN3
      IF (IFLAG .EQ. 0) RETURN
C
C     Should we evaluate functions or Jacobian?
C
      IF (IFLAG .EQ. 1) THEN
C
C       Evaluate functions.
C
        DO 10 I = 1,M
          TEMP = I
          FVEC(I) = TWO + TWO*TEMP - EXP(TEMP*X(1)) - EXP(TEMP*X(2))
   10   CONTINUE
      ELSE
C
C       Evaluate one row of Jacobian.
C
        IF (IFLAG .NE. 3) RETURN
        TEMP = NROW
        FJROW(1) = -TEMP*EXP(TEMP*X(1))
        FJROW(2) = -TEMP*EXP(TEMP*X(2))
      ENDIF
      RETURN
      END
