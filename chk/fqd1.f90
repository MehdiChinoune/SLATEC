!*==FQD1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FQD1
REAL FUNCTION FQD1(X)
  IMPLICIT NONE
  !*--FQD15
  !***BEGIN PROLOGUE  FQD1
  !***SUBSIDIARY
  !***PURPOSE  Function evaluator for QNC79 and GAUS8 quick checks.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FQD1-S, DFQD1-D)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***SEE ALSO  QG8TST, QN79QX
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   920229  DATE WRITTEN
  !***END PROLOGUE  FQD1
  !     .. Scalar Arguments ..
  REAL X
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !***FIRST EXECUTABLE STATEMENT  FQD1
  FQD1 = 0.0E0
  IF ( X>0.0E0 ) FQD1 = 1.0E0/SQRT(X)
END FUNCTION FQD1
