!** COSDG
REAL FUNCTION COSDG(X)
  !>
  !***
  !  Compute the cosine of an argument in degrees.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      SINGLE PRECISION (COSDG-S, DCOSDG-D)
  !***
  ! **Keywords:**  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
  !             TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! COSDG(X) evaluates the cosine for real X in degrees.
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

  ! JUNE 1977 EDITION.   W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
  REAL X
  INTEGER n
  REAL, PARAMETER :: raddeg = .017453292519943296E0
  !
  !* FIRST EXECUTABLE STATEMENT  COSDG
  COSDG = COS(raddeg*X)
  !
  IF ( MOD(X,90.)/=0. ) RETURN
  n = INT( ABS(X)/90.0 + 0.5 )
  n = MOD(n,2)
  IF ( n==0 ) COSDG = SIGN(1.0,COSDG)
  IF ( n==1 ) COSDG = 0.0
  !
END FUNCTION COSDG
