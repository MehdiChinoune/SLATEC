!DECK DCOSDG
REAL(8) FUNCTION DCOSDG(X)
  IMPLICIT NONE
  INTEGER n
  !***BEGIN PROLOGUE  DCOSDG
  !***PURPOSE  Compute the cosine of an argument in degrees.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      DOUBLE PRECISION (COSDG-S, DCOSDG-D)
  !***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
  !             TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DCOSDG(X) calculates the double precision trigonometric cosine
  ! for double precision argument X in units of degrees.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DCOSDG
  REAL(8) :: X, raddeg
  SAVE raddeg
  DATA raddeg/0.017453292519943295769236907684886D0/
  !***FIRST EXECUTABLE STATEMENT  DCOSDG
  DCOSDG = COS(raddeg*X)
  !
  IF ( MOD(X,90.D0)/=0.D0 ) RETURN
  n = INT( ABS(X)/90.D0 + 0.5D0 )
  n = MOD(n,2)
  IF ( n==0 ) DCOSDG = SIGN(1.0D0,DCOSDG)
  IF ( n==1 ) DCOSDG = 0.0D0
  !
END FUNCTION DCOSDG
