!DECK COSDG
REAL FUNCTION COSDG(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  COSDG
  !***PURPOSE  Compute the cosine of an argument in degrees.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      SINGLE PRECISION (COSDG-S, DCOSDG-D)
  !***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
  !             TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! COSDG(X) evaluates the cosine for real X in degrees.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  COSDG
  ! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
  REAL raddeg, X
  INTEGER n
  SAVE raddeg
  DATA raddeg/.017453292519943296E0/
  !
  !***FIRST EXECUTABLE STATEMENT  COSDG
  COSDG = COS(raddeg*X)
  !
  IF ( MOD(X,90.)/=0. ) RETURN
  n = INT( ABS(X)/90.0 + 0.5 )
  n = MOD(n,2)
  IF ( n==0 ) COSDG = SIGN(1.0,COSDG)
  IF ( n==1 ) COSDG = 0.0
  !
END FUNCTION COSDG
