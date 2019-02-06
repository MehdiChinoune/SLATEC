!*==F5S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F5S
      REAL FUNCTION F5S(X)
      IMPLICIT NONE
!*--F5S5
!***BEGIN PROLOGUE  F5S
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F5S
      REAL X
!***FIRST EXECUTABLE STATEMENT  F5S
      F5S = 0.0
      IF ( X/=0.0 ) F5S = 1.0/(X*SQRT(X))
      END FUNCTION F5S
