!*==DDRES1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DDRES1
      SUBROUTINE DDRES1(T,Y,Yprime,Delta,Ires,Rpar,Ipar)
      IMPLICIT NONE
!*--DDRES15
!***BEGIN PROLOGUE  DDRES1
!***SUBSIDIARY
!***PURPOSE  First residual evaluator for DDASQC.
!***LIBRARY   SLATEC (DASSL)
!***TYPE      DOUBLE PRECISION (SDRES1-S, DDRES1-D)
!***AUTHOR  PETZOLD, LINDA R., (LLNL)
!***SEE ALSO  DDASQC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   891013  DATE WRITTEN
!   901001  Converted prologue to 4.0 format and made all argument
!           declarations explicit.  (FNF)
!***END PROLOGUE  DDRES1
      INTEGER Ires , Ipar(*)
      DOUBLE PRECISION T , Y(*) , Yprime(*) , Delta(*) , Rpar(*)
!***FIRST EXECUTABLE STATEMENT  DDRES1
      Delta(1) = Yprime(1) + 10.0D0*Y(1)
      Delta(2) = Y(2) + Y(1) - 1.0D0
      END SUBROUTINE DDRES1
