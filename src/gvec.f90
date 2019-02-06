!*==GVEC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK GVEC
      SUBROUTINE GVEC(X,G)
      IMPLICIT NONE
!*--GVEC5
!*** Start of declarations inserted by SPAG
      REAL G , X
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  GVEC
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  GVEC
      DIMENSION G(*)
!***FIRST EXECUTABLE STATEMENT  GVEC
      G(1) = 0.0
      G(2) = 1.0 + COS(X)
      END SUBROUTINE GVEC
