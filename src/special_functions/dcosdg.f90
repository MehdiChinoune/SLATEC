!** DCOSDG
REAL(DP) ELEMENTAL FUNCTION DCOSDG(X)
  !> Compute the cosine of an argument in degrees.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      DOUBLE PRECISION (COSDG-S, DCOSDG-D)
  !***
  ! **Keywords:**  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DCOSDG(X) calculates the double precision trigonometric cosine
  ! for double precision argument X in units of degrees.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP), INTENT(IN) :: X
  INTEGER :: n
  REAL(DP), PARAMETER :: raddeg = 0.017453292519943295769236907684886_DP
  !* FIRST EXECUTABLE STATEMENT  DCOSDG
  DCOSDG = COS(raddeg*X)
  !
  IF( MOD(X,90._DP)/=0._DP ) RETURN
  n = INT( ABS(X)/90._DP + 0.5_DP )
  n = MOD(n,2)
  IF( n==0 ) DCOSDG = SIGN(1._DP,DCOSDG)
  IF( n==1 ) DCOSDG = 0._DP
  !
END FUNCTION DCOSDG