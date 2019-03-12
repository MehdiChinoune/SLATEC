!DECK EI
FUNCTION EI(X)
  IMPLICIT NONE
  REAL E1, EI, X
  !***BEGIN PROLOGUE  EI
  !***PURPOSE  Compute the exponential integral Ei(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C5
  !***TYPE      SINGLE PRECISION (EI-S, DEI-D)
  !***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! EI calculates the single precision exponential integral, Ei(X), for
  ! positive single precision argument X and the Cauchy principal value
  ! for negative X.  If principal values are used everywhere, then, for
  ! all X,
  !
  !    Ei(X) = -E1(-X)
  ! or
  !    E1(X) = -Ei(-X).
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  E1
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   891115  Modified prologue description.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  EI
  !***FIRST EXECUTABLE STATEMENT  EI
  EI = -E1(-X)
  !
END FUNCTION EI
