!*==DF1F.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF1F
      DOUBLE PRECISION FUNCTION DF1F(X)
      IMPLICIT NONE
!*--DF1F5
!***BEGIN PROLOGUE  DF1F
!***PURPOSE  Subsidiary to
!***LIBRARY   SLATEC
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DF1F
      DOUBLE PRECISION X , x1 , y
!***FIRST EXECUTABLE STATEMENT  DF1F
      x1 = X + 0.1D+01
      DF1F = 0.5D+01/x1/x1
      y = 0.5D+01/x1
      IF ( y>3.1415926535897932D0 ) DF1F = 0.0D0
      END FUNCTION DF1F
