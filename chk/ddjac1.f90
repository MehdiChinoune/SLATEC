!*==DDJAC1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DDJAC1
SUBROUTINE DDJAC1(T,Y,Yprime,Pd,Cj,Rpar,Ipar)
  IMPLICIT NONE
  !*--DDJAC15
  !***BEGIN PROLOGUE  DDJAC1
  !***SUBSIDIARY
  !***PURPOSE  First Jacobian evaluator for DDASQC.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      DOUBLE PRECISION (SDJAC1-S, DDJAC1-D)
  !***AUTHOR  PETZOLD, LINDA R., (LLNL)
  !***SEE ALSO  DDASQC
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   891013  DATE WRITTEN
  !   901001  Converted prologue to 4.0 format and made all argument
  !           declarations explicit.  (FNF)
  !***END PROLOGUE  DDJAC1
  INTEGER Ipar(*)
  DOUBLE PRECISION T , Y(*) , Yprime(*) , Pd(2,2) , Cj , Rpar(*)
  !***FIRST EXECUTABLE STATEMENT  DDJAC1
  Pd(1,1) = Cj + 10.0D0
  Pd(1,2) = 0.0D0
  Pd(2,1) = 1.0D0
  Pd(2,2) = 1.0D0
END SUBROUTINE DDJAC1
