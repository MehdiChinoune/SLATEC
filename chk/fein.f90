!*==FEIN.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK FEIN
      REAL FUNCTION FEIN(T)
      IMPLICIT NONE
!*--FEIN5
!***BEGIN PROLOGUE  FEIN
!***PURPOSE  Subsidiary to EG8CK.
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    FEINX
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FEIN
      COMMON /FEINX / X , A , FKM
      REAL X , A , FKM , T , aln
!***FIRST EXECUTABLE STATEMENT  FEIN
      aln = (FKM-T)*X - A*LOG(T)
      FEIN = EXP(aln)
      END FUNCTION FEIN
