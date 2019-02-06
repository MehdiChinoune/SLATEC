!*==SOSFNC.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SOSFNC
      REAL FUNCTION SOSFNC(X,K)
      IMPLICIT NONE
!*--SOSFNC5
!*** Start of declarations inserted by SPAG
      INTEGER K
      REAL X
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  SOSFNC
!***PURPOSE  Function evaluator for SOS quick check.
!***LIBRARY   SLATEC
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     FUNCTION WHICH EVALUATES THE FUNCTIONS, ONE AT A TIME,
!     FOR TEST PROGRAM USED IN QUICK CHECK OF SOS.
!
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890618  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  SOSFNC
      DIMENSION X(2)
!***FIRST EXECUTABLE STATEMENT  SOSFNC
      IF ( K==1 ) SOSFNC = 1.E0 - X(1)
      IF ( K==2 ) SOSFNC = 1.E1*(X(2)-X(1)**2)
      END FUNCTION SOSFNC
