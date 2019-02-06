*DECK FDTRUE
      SUBROUTINE FDTRUE (X, F, D)
C***BEGIN PROLOGUE  FDTRUE
C***SUBSIDIARY
C***PURPOSE  Compute exact function values for EVCHCK.
C***LIBRARY   SLATEC (PCHIP)
C***TYPE      SINGLE PRECISION (FDTRUE-S, DFDTRU-D)
C***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
C***AUTHOR  Fritsch, F. N., (LLNL)
C***DESCRIPTION
C
C        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
C
C                   F(X) = X*(X+1)*(X-2)
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890618  REVISION DATE from Version 3.2
C   890706  Cosmetic changes to prologue.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Revised prologue.  (FNF)
C   900316  Deleted variables ONE and TWO.  (FNF)
C   900321  Changed name of d.p. version from DFTRUE to DFDTRU.
C***END PROLOGUE  FDTRUE
      REAL  X, F, D
      DOUBLE PRECISION  FACT1, FACT2, XX
C
C***FIRST EXECUTABLE STATEMENT  FDTRUE
      XX = X
      FACT1 = XX + 1
      FACT2 = XX - 2
      F = XX * FACT1 * FACT2
      D = FACT1*FACT2 + XX*(FACT1 + FACT2)
C
      RETURN
C------------- LAST LINE OF FDTRUE FOLLOWS -----------------------------
      END
