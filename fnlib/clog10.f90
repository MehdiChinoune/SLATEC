!*==CLOG10.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK CLOG10
      COMPLEX FUNCTION CLOG10(Z)
      IMPLICIT NONE
!*--CLOG105
!*** Start of declarations inserted by SPAG
      REAL aloge
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CLOG10
!***PURPOSE  Compute the principal value of the complex base 10
!            logarithm.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      COMPLEX (CLOG10-C)
!***KEYWORDS  BASE TEN LOGARITHM, ELEMENTARY FUNCTIONS, FNLIB
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CLOG10(Z) calculates the principal value of the complex common
! or base 10 logarithm of Z for -PI .LT. arg(Z) .LE. +PI.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CLOG10
      COMPLEX Z
      SAVE aloge
      DATA aloge/0.43429448190325182765E0/
!***FIRST EXECUTABLE STATEMENT  CLOG10
      CLOG10 = aloge*LOG(Z)
!
      END FUNCTION CLOG10
