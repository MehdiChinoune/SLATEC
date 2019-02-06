*DECK SQFCN2
      SUBROUTINE SQFCN2 (N, X, FVEC, IFLAG)
C***BEGIN PROLOGUE  SQFCN2
C***PURPOSE  Evaluate function used in SNSQE.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SQFCN2-S, DQFCN2-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Subroutine which evaluates the function for test program
C   used in quick check of SNSQE.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  TYPE and declarations sections added.  (WRB)
C***END PROLOGUE  SQFCN2
C     .. Scalar Arguments ..
      INTEGER IFLAG, N
C     .. Array Arguments ..
      REAL FVEC(*), X(*)
C***FIRST EXECUTABLE STATEMENT  SQFCN2
      FVEC(1) = 1.0E0 - X(1)
      FVEC(2) = 10.0E0*(X(2)-X(1)**2)
      RETURN
      END
