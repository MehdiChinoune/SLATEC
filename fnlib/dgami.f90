!DECK DGAMI
REAL(8) FUNCTION DGAMI(A,X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DGAMI
  !***PURPOSE  Evaluate the incomplete Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      DOUBLE PRECISION (GAMI-S, DGAMI-D)
  !***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate the incomplete gamma function defined by
  !
  ! DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
  !
  ! DGAMI is evaluated for positive values of A and non-negative values
  ! of X.  A slight deterioration of 2 or 3 digits accuracy will occur
  ! when DGAMI is very large or very small, because logarithmic variables
  ! are used.  The function and both arguments are double precision.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DGAMIT, DLNGAM, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DGAMI
  REAL(8) :: A, X, factor, DLNGAM, DGAMIT
  !***FIRST EXECUTABLE STATEMENT  DGAMI
  IF ( A<=0.D0 ) CALL XERMSG('SLATEC','DGAMI','A MUST BE GT ZERO',1,2)
  IF ( X<0.D0 ) CALL XERMSG('SLATEC','DGAMI','X MUST BE GE ZERO',2,2)
  !
  DGAMI = 0.D0
  IF ( X==0.0D0 ) RETURN
  !
  ! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
  factor = EXP(DLNGAM(A)+A*LOG(X))
  !
  DGAMI = factor*DGAMIT(A,X)
  !
END FUNCTION DGAMI
