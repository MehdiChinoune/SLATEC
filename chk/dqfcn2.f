*DECK DQFCN2
      SUBROUTINE DQFCN2 (N, X, FVEC, IFLAG)
C***BEGIN PROLOGUE  DQFCN2
C***PURPOSE  Evaluate function used in DNSQE.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SQFCN2-S, DQFCN2-D)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Subroutine which evaluates the function for test program
C   used in quick check of DNSQE.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   930214  TYPE and declarations sections added.  (WRB)
C***END PROLOGUE  DQFCN2
C     .. Scalar Arguments ..
      INTEGER IFLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION FVEC(*), X(*)
C***FIRST EXECUTABLE STATEMENT  DQFCN2
      FVEC(1) = 1.0D0 - X(1)
      FVEC(2) = 10.0D0*(X(2)-X(1)**2)
      RETURN
      END
