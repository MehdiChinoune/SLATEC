!*==DFQD1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFQD1
REAL(8) FUNCTION DFQD1(X)
  IMPLICIT NONE
  !*--DFQD15
  !***BEGIN PROLOGUE  DFQD1
  !***SUBSIDIARY
  !***PURPOSE  Function evaluator for DQNC79 and DGAUS8 quick checks.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FQD1-S, DFQD1-D)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***SEE ALSO  DQG8TS, DQN79Q
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   920229  DATE WRITTEN
  !***END PROLOGUE  DFQD1
  !     .. Scalar Arguments ..
  REAL(8) :: X
  !     .. Intrinsic Functions ..
  INTRINSIC SQRT
  !***FIRST EXECUTABLE STATEMENT  DFQD1
  DFQD1 = 0.0D0
  IF ( X>0.0D0 ) DFQD1 = 1.0D0/SQRT(X)
END FUNCTION DFQD1
