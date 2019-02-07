!*==GAMI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK GAMI
FUNCTION GAMI(A,X)
  IMPLICIT NONE
  !*--GAMI5
  !*** Start of declarations inserted by SPAG
  REAL A , ALNGAM , factor , GAMI , GAMIT , X
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  GAMI
  !***PURPOSE  Evaluate the incomplete Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      SINGLE PRECISION (GAMI-S, DGAMI-D)
  !***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate the incomplete gamma function defined by
  !
  ! GAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0) .
  !
  ! GAMI is evaluated for positive values of A and non-negative values
  ! of X.  A slight deterioration of 2 or 3 digits accuracy will occur
  ! when GAMI is very large or very small, because logarithmic variables
  ! are used.  GAMI, A, and X are single precision.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ALNGAM, GAMIT, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  GAMI
  !***FIRST EXECUTABLE STATEMENT  GAMI
  IF ( A<=0.0 ) CALL XERMSG('SLATEC','GAMI','A MUST BE GT ZERO',1,2)
  IF ( X<0.0 ) CALL XERMSG('SLATEC','GAMI','X MUST BE GE ZERO',2,2)
  !
  GAMI = 0.0
  IF ( X==0.0 ) RETURN
  !
  ! THE ONLY ERROR POSSIBLE IN THE EXPRESSION BELOW IS A FATAL OVERFLOW.
  factor = EXP(ALNGAM(A)+A*LOG(X))
  !
  GAMI = factor*GAMIT(A,X)
  !
END FUNCTION GAMI
