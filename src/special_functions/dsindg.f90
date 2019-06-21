!** DSINDG
REAL(DP) FUNCTION DSINDG(X)
  !> Compute the sine of an argument in degrees.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      DOUBLE PRECISION (SINDG-S, DSINDG-D)
  !***
  ! **Keywords:**  DEGREES, ELEMENTARY FUNCTIONS, FNLIB, SINE, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DSINDG(X) calculates the double precision sine for double
  ! precision argument X where X is in degrees.
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

  INTEGER :: n
  REAL(DP) :: X
  REAL(DP), PARAMETER :: raddeg = 0.017453292519943295769236907684886_DP
  !* FIRST EXECUTABLE STATEMENT  DSINDG
  DSINDG = SIN(raddeg*X)
  !
  IF( MOD(X,90._DP)/=0._DP ) RETURN
  n = INT( ABS(X)/90._DP + 0.5_DP )
  n = MOD(n,2)
  IF( n==0 ) DSINDG = 0._DP
  IF( n==1 ) DSINDG = SIGN(1._DP,DSINDG)
  !
END FUNCTION DSINDG
