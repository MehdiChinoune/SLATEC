!** R9ATN1
REAL FUNCTION R9ATN1(X)
  !>
  !  Evaluate ATAN(X) from first order relative accuracy so that
  !            ATAN(X) = X + X**3*R9ATN1(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      SINGLE PRECISION (R9ATN1-S, D9ATN1-D)
  !***
  ! **Keywords:**  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB,
  !             TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  ATAN(X)  from first order, that is, evaluate
  ! (ATAN(X)-X)/X**3  with relative error accuracy so that
  !        ATAN(X) = X + X**3*R9ATN1(X).
  !
  ! Series for ATN1       on the interval  0.          to  1.00000D+00
  !                                        with weighted error   2.21E-17
  !                                         log weighted error  16.66
  !                               significant figures required  15.44
  !                                    decimal places required  17.32
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL eps, X, y
  INTEGER, SAVE :: ntatn1
  REAL, SAVE :: xsml, xbig, xmax
  REAL, PARAMETER :: atn1cs(21) = [ -.03283997535355202E0, .05833432343172412E0, &
    -.00740036969671964E0, .00100978419933728E0,-.00014397871635652E0, &
    .00002114512648992E0, -.00000317232107425E0, .00000048366203654E0, &
    -.00000007467746546E0, .00000001164800896E0,-.00000000183208837E0, &
    .00000000029019082E0, -.00000000004623885E0, .00000000000740552E0, &
    -.00000000000119135E0, .00000000000019240E0,-.00000000000003118E0, &
    .00000000000000506E0, -.00000000000000082E0, .00000000000000013E0, &
    -.00000000000000002E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  R9ATN1
  IF ( first ) THEN
    eps = R1MACH(3)
    ntatn1 = INITS(atn1cs,21,0.1*eps)
    !
    xsml = SQRT(0.1*eps)
    xbig = 1.571/SQRT(eps)
    xmax = 1.571/eps
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF ( y>1.0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','R9ATN1',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',2,2)
    IF ( y>xbig ) CALL XERMSG('SLATEC','R9ATN1',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG',1,1)
    !
    R9ATN1 = (ATAN(X)-X)/X**3
    RETURN
  END IF
  !
  IF ( y<=xsml ) R9ATN1 = -1.0/3.0
  IF ( y<=xsml ) RETURN
  !
  R9ATN1 = -0.25 + CSEVL(2.0*y*y-1.,atn1cs,ntatn1)
  RETURN
END FUNCTION R9ATN1
