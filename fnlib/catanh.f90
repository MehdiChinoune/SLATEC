!*==CATANH.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK CATANH
      COMPLEX FUNCTION CATANH(Z)
      IMPLICIT NONE
!*--CATANH5
!***BEGIN PROLOGUE  CATANH
!***PURPOSE  Compute the arc hyperbolic tangent.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      COMPLEX (ATANH-S, DATANH-D, CATANH-C)
!***KEYWORDS  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS,
!             FNLIB, INVERSE HYPERBOLIC TANGENT
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CATANH(Z) calculates the complex arc hyperbolic tangent of Z.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CATAN
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CATANH
      COMPLEX Z , ci , CATAN
      SAVE ci
      DATA ci/(0.,1.)/
!***FIRST EXECUTABLE STATEMENT  CATANH
      CATANH = -ci*CATAN(ci*Z)
!
      END FUNCTION CATANH
