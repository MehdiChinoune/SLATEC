!*==F4P.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F4P
      REAL FUNCTION F4P(X)
      IMPLICIT NONE
!*--F4P5
!***BEGIN PROLOGUE  F4P
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F4P
      REAL X
!***FIRST EXECUTABLE STATEMENT  F4P
      F4P = 0.0
      IF ( X>0.0 ) F4P = 1.0/(X*SQRT(X))
      END FUNCTION F4P
