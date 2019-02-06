!*==F1WS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1WS
      REAL FUNCTION F1WS(X)
      IMPLICIT NONE
!*--F1WS5
!***BEGIN PROLOGUE  F1WS
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F1WS
      REAL X
!***FIRST EXECUTABLE STATEMENT  F1WS
      F1WS = ABS(X-0.33E+00)**(-0.999E+00)
      END FUNCTION F1WS
