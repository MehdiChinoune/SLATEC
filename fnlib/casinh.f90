!*==CASINH.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK CASINH
      COMPLEX FUNCTION CASINH(Z)
      IMPLICIT NONE
!*--CASINH5
!***BEGIN PROLOGUE  CASINH
!***PURPOSE  Compute the arc hyperbolic sine.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      COMPLEX (ASINH-S, DASINH-D, CASINH-C)
!***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB,
!             INVERSE HYPERBOLIC SINE
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CASINH(Z) calculates the complex arc hyperbolic sine of Z.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CASIN
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CASINH
      COMPLEX Z , ci , CASIN
      SAVE ci
      DATA ci/(0.,1.)/
!***FIRST EXECUTABLE STATEMENT  CASINH
      CASINH = -ci*CASIN(ci*Z)
!
      END FUNCTION CASINH
