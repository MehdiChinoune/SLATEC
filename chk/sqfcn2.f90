!DECK SQFCN2
SUBROUTINE SQFCN2(N,X,Fvec,Iflag)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SQFCN2
  !***PURPOSE  Evaluate function used in SNSQE.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (SQFCN2-S, DQFCN2-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine which evaluates the function for test program
  !   used in quick check of SNSQE.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added.  (WRB)
  !***END PROLOGUE  SQFCN2
  !     .. Scalar Arguments ..
  INTEGER Iflag, N
  !     .. Array Arguments ..
  REAL Fvec(*), X(*)
  !***FIRST EXECUTABLE STATEMENT  SQFCN2
  Fvec(1) = 1.0E0 - X(1)
  Fvec(2) = 10.0E0*(X(2)-X(1)**2)
END SUBROUTINE SQFCN2
