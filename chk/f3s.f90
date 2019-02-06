!*==F3S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F3S
      REAL FUNCTION F3S(X)
      IMPLICIT NONE
!*--F3S5
!***BEGIN PROLOGUE  F3S
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F3S
      REAL X
!***FIRST EXECUTABLE STATEMENT  F3S
      F3S = 0.1E+01
      IF ( X>3.1415926535897932 ) F3S = 0.0
      END FUNCTION F3S
