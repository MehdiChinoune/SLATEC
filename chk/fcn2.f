*DECK FCN2
      SUBROUTINE FCN2 (IFLAG, M, N, X, FVEC, FJAC, LDFJAC)
C***BEGIN PROLOGUE  FCN2
C***PURPOSE  Subsidiary to SNLS1Q.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FCN2-S, DFCN2-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Subroutine to evaluate function and full Jacobian for test
C   problem in quick check of SNLS1E.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  TYPE and declarations sections added and code polished.
C           (WRB)
C***END PROLOGUE  FCN2
C     .. Scalar Arguments ..
      INTEGER IFLAG, LDFJAC, M, N
C     .. Array Arguments ..
      REAL FJAC(LDFJAC,*), FVEC(*), X(*)
C     .. Local Scalars ..
      REAL TEMP, TWO
      INTEGER I
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     .. Data statements ..
      DATA TWO /2.0E0/
C***FIRST EXECUTABLE STATEMENT  FCN2
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
C       Evaluate Jacobian.
C
        IF (IFLAG .NE. 2) RETURN
        DO 20 I = 1,M
          TEMP = I
          FJAC(I,1) = -TEMP*EXP(TEMP*X(1))
          FJAC(I,2) = -TEMP*EXP(TEMP*X(2))
   20   CONTINUE
      ENDIF
      RETURN
      END
