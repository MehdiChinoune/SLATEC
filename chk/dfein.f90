!*==DFEIN.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DFEIN
REAL(8) FUNCTION DFEIN(T)
  IMPLICIT NONE
  !*--DFEIN5
  !***BEGIN PROLOGUE  DFEIN
  !***PURPOSE  Subsidiary to DEG8CK.
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    DFEINX
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DFEIN
  COMMON /DFEINX/ X , A , FKM
  REAL(8) :: X , A , FKM , T , aln
  !***FIRST EXECUTABLE STATEMENT  DFEIN
  aln = (FKM-T)*X - A*LOG(T)
  DFEIN = EXP(aln)
END FUNCTION DFEIN
