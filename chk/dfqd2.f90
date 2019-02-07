!*==DFQD2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFQD2
REAL(8) FUNCTION DFQD2(X)
  IMPLICIT NONE
  !*--DFQD25
  !***BEGIN PROLOGUE  DFQD2
  !***SUBSIDIARY
  !***PURPOSE  Function evaluator for DQNC79 and DGAUS8 quick checks.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (FQD2-S, DFQD2-D)
  !***AUTHOR  Boland, W. Robert, (LANL)
  !***SEE ALSO  DQG8TS, DQN79Q
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   920229  DATE WRITTEN
  !***END PROLOGUE  DFQD2
  !     .. Scalar Arguments ..
  REAL(8) :: X
  !     .. Intrinsic Functions ..
  INTRINSIC COS, EXP
  !***FIRST EXECUTABLE STATEMENT  DFQD2
  DFQD2 = EXP(X)*COS(10.0D0*X)
END FUNCTION DFQD2
