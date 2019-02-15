!DECK DQFCN2
SUBROUTINE DQFCN2(N,X,Fvec,Iflag)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DQFCN2
  !***PURPOSE  Evaluate function used in DNSQE.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SQFCN2-S, DQFCN2-D)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !   Subroutine which evaluates the function for test program
  !   used in quick check of DNSQE.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   930214  TYPE and declarations sections added.  (WRB)
  !***END PROLOGUE  DQFCN2
  !     .. Scalar Arguments ..
  INTEGER Iflag, N
  !     .. Array Arguments ..
  REAL(8) :: Fvec(*), X(*)
  !***FIRST EXECUTABLE STATEMENT  DQFCN2
  Fvec(1) = 1.0D0 - X(1)
  Fvec(2) = 10.0D0*(X(2)-X(1)**2)
END SUBROUTINE DQFCN2
