*DECK FDEQC
      SUBROUTINE FDEQC (T, U, UPRIME, RPAR, IPAR)
C***BEGIN PROLOGUE  FDEQC
C***SUBSIDIARY
C***PURPOSE  Derivative evaluator for DEPAC quick checks.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FDEQC-S, DFDEQC-D)
C***AUTHOR  Chow, Jeff, (LANL)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900415  Name changed from F to FDEQC.  (WRB)
C***END PROLOGUE  FDEQC
C
C     Declare arguments.
C
      INTEGER IPAR(*)
      REAL RPAR(*), T, U(*), UPRIME(*)
C
C     Declare local variables.
C
      REAL R, RSQ, R3
C***FIRST EXECUTABLE STATEMENT  FDEQC
      RSQ = U(1)*U(1) + U(2)*U(2)
      R = SQRT(RSQ)
      R3 = RSQ*R
      UPRIME(1) = U(3)
      UPRIME(2) = U(4)
      UPRIME(3) = -(U(1)/R3)
      UPRIME(4) = -(U(2)/R3)
      RETURN
      END
