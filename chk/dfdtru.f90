!*==DFDTRU.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFDTRU
SUBROUTINE DFDTRU(X,F,D)
  IMPLICIT NONE
  !*--DFDTRU5
  !***BEGIN PROLOGUE  DFDTRU
  !***SUBSIDIARY
  !***PURPOSE  Compute exact function values for DEVCHK.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (FDTRUE-S, DFDTRU-D)
  !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
  !***AUTHOR  Fritsch, F. N., (LLNL)
  !***DESCRIPTION
  !
  !        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
  !
  !                   F(X) = X*(X+1)*(X-2)
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   890706  Cosmetic changes to prologue.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  Revised prologue.  (FNF)
  !   900316  Deleted variables ONE and TWO.  (FNF)
  !   900321  Changed name of d.p. version from DFTRUE to DFDTRU.
  !***END PROLOGUE  DFDTRU
  DOUBLE PRECISION X , F , D
  DOUBLE PRECISION fact1 , fact2 , xx
  !
  !***FIRST EXECUTABLE STATEMENT  DFDTRU
  xx = X
  fact1 = xx + 1
  fact2 = xx - 2
  F = xx*fact1*fact2
  D = fact1*fact2 + xx*(fact1+fact2)
  !
  !------------- LAST LINE OF DFDTRU FOLLOWS -----------------------------
END SUBROUTINE DFDTRU
