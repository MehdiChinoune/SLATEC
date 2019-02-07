!*==SINDG.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK SINDG
FUNCTION SINDG(X)
  IMPLICIT NONE
  !*--SINDG5
  !*** Start of declarations inserted by SPAG
  INTEGER n
  REAL raddeg, SINDG, X
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SINDG
  !***PURPOSE  Compute the sine of an argument in degrees.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      SINGLE PRECISION (SINDG-S, DSINDG-D)
  !***KEYWORDS  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! SINDG(X) evaluates the single precision sine of X where
  ! X is in degrees.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  SINDG
  ! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
  SAVE raddeg
  DATA raddeg/.017453292519943296E0/
  !
  !***FIRST EXECUTABLE STATEMENT  SINDG
  SINDG = SIN(raddeg*X)
  !
  IF ( MOD(X,90.)/=0. ) RETURN
  n = ABS(X)/90.0 + 0.5
  n = MOD(n,2)
  IF ( n==0 ) SINDG = 0.
  IF ( n==1 ) SINDG = SIGN(1.0,SINDG)
  !
END FUNCTION SINDG
