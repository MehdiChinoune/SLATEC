!*==F3G.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F3G
      REAL FUNCTION F3G(X)
      IMPLICIT NONE
!*--F3G5
!***BEGIN PROLOGUE  F3G
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F3G
      REAL X
!***FIRST EXECUTABLE STATEMENT  F3G
      F3G = ABS(X-0.33E+00)**(-0.9E+00)
      END FUNCTION F3G
