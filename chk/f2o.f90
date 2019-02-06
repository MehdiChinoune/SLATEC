!*==F2O.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F2O
      REAL FUNCTION F2O(X)
      IMPLICIT NONE
!*--F2O5
!***BEGIN PROLOGUE  F2O
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F2O
      REAL X
!***FIRST EXECUTABLE STATEMENT  F2O
      F2O = 0.0E+00
      IF ( X/=0.0E+00 ) F2O = 1.0/(X*X*SQRT(X))
      END FUNCTION F2O
