!*==F1O.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1O
      REAL FUNCTION F1O(X)
      IMPLICIT NONE
!*--F1O5
!***BEGIN PROLOGUE  F1O
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F1O
      REAL X
!***FIRST EXECUTABLE STATEMENT  F1O
      F1O = 1.0
      IF ( X>3.1415926535897932 ) F1O = 0.0
      END FUNCTION F1O
