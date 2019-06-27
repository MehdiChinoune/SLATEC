!** DEI
REAL(DP) ELEMENTAL FUNCTION DEI(X)
  !> Compute the exponential integral Ei(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      DOUBLE PRECISION (EI-S, DEI-D)
  !***
  ! **Keywords:**  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DE1

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   891115  Modified prologue description.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP), INTENT(IN) :: X
  !* FIRST EXECUTABLE STATEMENT  DEI
  DEI = -DE1(-X)
  !
END FUNCTION DEI