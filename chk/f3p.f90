!*==F3P.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F3P
      REAL FUNCTION F3P(X)
      IMPLICIT NONE
!*--F3P5
!***BEGIN PROLOGUE  F3P
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  F3P
      REAL X
!***FIRST EXECUTABLE STATEMENT  F3P
      F3P = 1.0
      IF ( X>3.1415926535897932 ) F3P = 0.0
      END FUNCTION F3P
