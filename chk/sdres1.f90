!*==SDRES1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SDRES1
SUBROUTINE SDRES1(T,Y,Yprime,Delta,Ires,Rpar,Ipar)
  IMPLICIT NONE
  !*--SDRES15
  !***BEGIN PROLOGUE  SDRES1
  !***SUBSIDIARY
  !***PURPOSE  First residual evaluator for SDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDRES1-S, DDRES1-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  SDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !***END PROLOGUE  SDRES1
  INTEGER Ires , Ipar(*)
  REAL T , Y(*) , Yprime(*) , Delta(*) , Rpar(*)
  !***FIRST EXECUTABLE STATEMENT  SDRES1
  Delta(1) = Yprime(1) + 10.0E0*Y(1)
  Delta(2) = Y(2) + Y(1) - 1.0E0
END SUBROUTINE SDRES1
