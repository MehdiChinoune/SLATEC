!** R9LGIT
REAL FUNCTION R9LGIT(A,X,Algap1)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithm of Tricomi's incomplete Gamma
  !            function with Perron's continued fraction for large X and
  !            A .GE. X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (R9LGIT-S, D9LGIT-D)
  !***
  ! **Keywords:**  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
  !             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the log of Tricomi's incomplete gamma function with Perron's
  ! continued fraction for large X and for A .GE. X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  
  REAL A, a1x, Algap1, ax, eps, fk, hstar, p, r, R1MACH, s, sqeps, t, X
  INTEGER k
  SAVE eps, sqeps
  DATA eps, sqeps/2*0.0/
  !* FIRST EXECUTABLE STATEMENT  R9LGIT
  IF ( eps==0.0 ) eps = 0.5*R1MACH(3)
  IF ( sqeps==0.0 ) sqeps = SQRT(R1MACH(4))
  !
  IF ( X<=0.0.OR.A<X ) CALL XERMSG('SLATEC','R9LGIT',&
    'X SHOULD BE GT 0.0 AND LE A',2,2)
  !
  ax = A + X
  a1x = ax + 1.0
  r = 0.0
  p = 1.0
  s = p
  DO k = 1, 200
    fk = k
    t = (A+fk)*X*(1.0+r)
    r = t/((ax+fk)*(a1x+fk)-t)
    p = r*p
    s = s + p
    IF ( ABS(p)<eps*s ) GOTO 100
  ENDDO
  CALL XERMSG('SLATEC','R9LGIT',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',3,2)
  !
  100  hstar = 1.0 - X*s/a1x
  IF ( hstar<sqeps ) CALL XERMSG('SLATEC','R9LGIT',&
    'RESULT LESS THAN HALF PRECISION',1,1)
  !
  R9LGIT = -X - Algap1 - LOG(hstar)
  !
END FUNCTION R9LGIT
