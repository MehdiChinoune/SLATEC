!*==DEI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DEI
      DOUBLE PRECISION FUNCTION DEI(X)
      IMPLICIT NONE
!*--DEI5
!***BEGIN PROLOGUE  DEI
!***PURPOSE  Compute the exponential integral Ei(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C5
!***TYPE      DOUBLE PRECISION (EI-S, DEI-D)
!***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DEI calculates the double precision exponential integral, Ei(X), for
! positive double precision argument X and the Cauchy principal value
! for negative X.  If principal values are used everywhere, then, for
! all X,
!
!    Ei(X) = -E1(-X)
! or
!    E1(X) = -Ei(-X).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DE1
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   891115  Modified prologue description.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DEI
      DOUBLE PRECISION X , DE1
!***FIRST EXECUTABLE STATEMENT  DEI
      DEI = -DE1(-X)
!
      END FUNCTION DEI
