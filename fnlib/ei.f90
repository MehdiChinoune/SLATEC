!*==EI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK EI
      FUNCTION EI(X)
      IMPLICIT NONE
!*--EI5
!*** Start of declarations inserted by SPAG
      REAL E1 , EI , X
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  EI
!***PURPOSE  Compute the exponential integral Ei(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      SINGLE PRECISION (EI-S, DEI-D)
!***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! EI calculates the single precision exponential integral, Ei(X), for
! positive single precision argument X and the Cauchy principal value
! for negative X.  If principal values are used everywhere, then, for
! all X,
!
!    Ei(X) = -E1(-X)
! or
!    E1(X) = -Ei(-X).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  E1
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   891115  Modified prologue description.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  EI
!***FIRST EXECUTABLE STATEMENT  EI
      EI = -E1(-X)
!
      END FUNCTION EI
