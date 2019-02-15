!DECK FQD2
REAL FUNCTION FQD2(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  FQD2
  !***SUBSIDIARY
  !***PURPOSE  Function evaluator for QNC79 and GAUS8 quick checks.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FQD2-S, DFQD2-D)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***SEE ALSO  QG8TST, QN79QX
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   920229  DATE WRITTEN
  !***END PROLOGUE  FQD2
  !     .. Scalar Arguments ..
  REAL X
  !     .. Intrinsic Functions ..
  INTRINSIC COS, EXP
  !***FIRST EXECUTABLE STATEMENT  FQD2
  FQD2 = EXP(X)*COS(10.0E0*X)
END FUNCTION FQD2
