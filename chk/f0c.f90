!*==F0C.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F0C
      REAL FUNCTION F0C(X)
      IMPLICIT NONE
!*--F0C5
!***BEGIN PROLOGUE  F0C
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F0C
      REAL X
!***FIRST EXECUTABLE STATEMENT  F0C
      F0C = 1.E0/(X*X+1.E-4)
      END FUNCTION F0C
