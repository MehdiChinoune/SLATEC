!*==ALGAMS.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK ALGAMS
      SUBROUTINE ALGAMS(X,Algam,Sgngam)
      IMPLICIT NONE
!*--ALGAMS5
!*** Start of declarations inserted by SPAG
      REAL Algam , ALNGAM , Sgngam , X
      INTEGER int
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  ALGAMS
!***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!            function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      SINGLE PRECISION (ALGAMS-S, DLGAMS-D)
!***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
!             FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluates the logarithm of the absolute value of the gamma
! function.
!     X           - input argument
!     ALGAM       - result
!     SGNGAM      - is set to the sign of GAMMA(X) and will
!                   be returned at +1.0 or -1.0.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNGAM
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ALGAMS
!***FIRST EXECUTABLE STATEMENT  ALGAMS
      Algam = ALNGAM(X)
      Sgngam = 1.0
      IF ( X>0.0 ) RETURN
!
      int = MOD(-AINT(X),2.0) + 0.1
      IF ( int==0 ) Sgngam = -1.0
!
      END SUBROUTINE ALGAMS
