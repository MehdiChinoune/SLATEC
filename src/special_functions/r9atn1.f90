!** R9ATN1
REAL(SP) ELEMENTAL FUNCTION R9ATN1(X)
  !> Evaluate ATAN(X) from first order relative accuracy so that
  !            ATAN(X) = X + X**3*R9ATN1(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      SINGLE PRECISION (R9ATN1-S, D9ATN1-D)
  !***
  ! **Keywords:**  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB, TRIGONOMETRIC
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
  USE service, ONLY : eps_2_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntatn1 = 10
  REAL(SP), PARAMETER :: eps = eps_2_sp, xsml = SQRT(0.1_SP*eps), &
    xbig = 1.571_SP/SQRT(eps), xmax = 1.571_SP/eps
  REAL(SP), PARAMETER :: atn1cs(21) = [ -.03283997535355202_SP, .05833432343172412_SP, &
    -.00740036969671964_SP, .00100978419933728_SP,-.00014397871635652_SP, &
    .00002114512648992_SP, -.00000317232107425_SP, .00000048366203654_SP, &
    -.00000007467746546_SP, .00000001164800896_SP,-.00000000183208837_SP, &
    .00000000029019082_SP, -.00000000004623885_SP, .00000000000740552_SP, &
    -.00000000000119135_SP, .00000000000019240_SP,-.00000000000003118_SP, &
    .00000000000000506_SP, -.00000000000000082_SP, .00000000000000013_SP, &
    -.00000000000000002_SP ]
  !* FIRST EXECUTABLE STATEMENT  R9ATN1
  ! ntatn1 = INITS(atn1cs,0.1_SP*eps)
  !
  y = ABS(X)
  IF( y>xmax ) THEN
    ERROR STOP 'R9ATN1 : NO PRECISION IN ANSWER BECAUSE X IS TOO BIG'
  ELSEIF( y>1._SP ) THEN
    ! IF( y>xbig ) 'R9ATN1 : ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'
    R9ATN1 = (ATAN(X)-X)/X**3
  ELSEIF( y>xsml ) THEN
    R9ATN1 = -0.25_SP + CSEVL(2._SP*y*y-1._SP,atn1cs(1:ntatn1))
  ELSE
    R9ATN1 = -1._SP/3._SP
  END IF

  RETURN
END FUNCTION R9ATN1