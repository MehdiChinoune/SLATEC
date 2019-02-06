!*==DLGAMS.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DLGAMS
      SUBROUTINE DLGAMS(X,Dlgam,Sgngam)
      IMPLICIT NONE
!*--DLGAMS5
!*** Start of declarations inserted by SPAG
      INTEGER int
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  DLGAMS
!***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!            function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (ALGAMS-S, DLGAMS-D)
!***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
!             FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLGAMS(X,DLGAM,SGNGAM) calculates the double precision natural
! logarithm of the absolute value of the Gamma function for
! double precision argument X and stores the result in double
! precision argument DLGAM.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DLNGAM
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DLGAMS
      DOUBLE PRECISION X , Dlgam , Sgngam , DLNGAM
!***FIRST EXECUTABLE STATEMENT  DLGAMS
      Dlgam = DLNGAM(X)
      Sgngam = 1.0D0
      IF ( X>0.D0 ) RETURN
!
      int = MOD(-AINT(X),2.0D0) + 0.1D0
      IF ( int==0 ) Sgngam = -1.0D0
!
      END SUBROUTINE DLGAMS
