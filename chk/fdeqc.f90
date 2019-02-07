!*==FDEQC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FDEQC
SUBROUTINE FDEQC(T,U,Uprime,Rpar,Ipar)
  IMPLICIT NONE
  !*--FDEQC5
  !***BEGIN PROLOGUE  FDEQC
  !***SUBSIDIARY
  !***PURPOSE  Derivative evaluator for DEPAC quick checks.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (FDEQC-S, DFDEQC-D)
  !***AUTHOR  Chow, Jeff, (LANL)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810801  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900415  Name changed from F to FDEQC.  (WRB)
  !***END PROLOGUE  FDEQC
  !
  !     Declare arguments.
  !
  INTEGER Ipar(*)
  REAL Rpar(*) , T , U(*) , Uprime(*)
  !
  !     Declare local variables.
  !
  REAL r , rsq , r3
  !***FIRST EXECUTABLE STATEMENT  FDEQC
  rsq = U(1)*U(1) + U(2)*U(2)
  r = SQRT(rsq)
  r3 = rsq*r
  Uprime(1) = U(3)
  Uprime(2) = U(4)
  Uprime(3) = -(U(1)/r3)
  Uprime(4) = -(U(2)/r3)
END SUBROUTINE FDEQC
