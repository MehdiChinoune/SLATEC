!*==SDJAC1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SDJAC1
SUBROUTINE SDJAC1(T,Y,Yprime,Pd,Cj,Rpar,Ipar)
  IMPLICIT NONE
  !*--SDJAC15
  !***BEGIN PROLOGUE  SDJAC1
  !***SUBSIDIARY
  !***PURPOSE  First Jacobian evaluator for SDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDJAC1-S, DDJAC1-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  SDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !***END PROLOGUE  SDJAC1
  INTEGER Ipar(*)
  REAL T, Y(*), Yprime(*), Pd(2,2), Cj, Rpar(*)
  !***FIRST EXECUTABLE STATEMENT  SDJAC1
  Pd(1,1) = Cj + 10.0E0
  Pd(1,2) = 0.0E0
  Pd(2,1) = 1.0E0
  Pd(2,2) = 1.0E0
END SUBROUTINE SDJAC1
