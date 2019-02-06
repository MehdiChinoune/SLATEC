*DECK DFCN1
      SUBROUTINE DFCN1 (IFLAG, M, N, X, FVEC, FJAC, LDFJAC)
C***BEGIN PROLOGUE  DFCN1
C***PURPOSE  Subsidiary to DNLS1Q.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (FCN1-S, DFCN1-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Subroutine which evaluates the function for test program
C   used in quick check of DNLS1E.
C
C   Numerical approximation of Jacobian is used.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  TYPE and declarations sections added.  (WRB)
C***END PROLOGUE  DFCN1
C     .. Scalar Arguments ..
      DOUBLE PRECISION FJAC
      INTEGER IFLAG, LDFJAC, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION FVEC(*), X(*)
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP, TWO
      INTEGER I
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     .. Data statements ..
      DATA TWO /2.0D0/
C***FIRST EXECUTABLE STATEMENT  DFCN1
      IF (IFLAG .NE. 1) RETURN
      DO 10 I = 1,M
        TEMP = I
        FVEC(I) = TWO + TWO*TEMP - EXP(TEMP*X(1)) - EXP(TEMP*X(2))
   10 CONTINUE
      RETURN
      END
