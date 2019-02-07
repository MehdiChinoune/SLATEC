!*==R9ATN1.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK R9ATN1
FUNCTION R9ATN1(X)
  IMPLICIT NONE
  !*--R9ATN15
  !*** Start of declarations inserted by SPAG
  REAL atn1cs , CSEVL , eps , R1MACH , R9ATN1 , X , xbig , xmax , xsml , y
  INTEGER INITS , ntatn1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  R9ATN1
  !***SUBSIDIARY
  !***PURPOSE  Evaluate ATAN(X) from first order relative accuracy so that
  !            ATAN(X) = X + X**3*R9ATN1(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      SINGLE PRECISION (R9ATN1-S, D9ATN1-D)
  !***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB,
  !             TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  R9ATN1
  DIMENSION atn1cs(21)
  LOGICAL first
  SAVE atn1cs , ntatn1 , xsml , xbig , xmax , first
  DATA atn1cs(1)/ - .03283997535355202E0/
  DATA atn1cs(2)/.05833432343172412E0/
  DATA atn1cs(3)/ - .00740036969671964E0/
  DATA atn1cs(4)/.00100978419933728E0/
  DATA atn1cs(5)/ - .00014397871635652E0/
  DATA atn1cs(6)/.00002114512648992E0/
  DATA atn1cs(7)/ - .00000317232107425E0/
  DATA atn1cs(8)/.00000048366203654E0/
  DATA atn1cs(9)/ - .00000007467746546E0/
  DATA atn1cs(10)/.00000001164800896E0/
  DATA atn1cs(11)/ - .00000000183208837E0/
  DATA atn1cs(12)/.00000000029019082E0/
  DATA atn1cs(13)/ - .00000000004623885E0/
  DATA atn1cs(14)/.00000000000740552E0/
  DATA atn1cs(15)/ - .00000000000119135E0/
  DATA atn1cs(16)/.00000000000019240E0/
  DATA atn1cs(17)/ - .00000000000003118E0/
  DATA atn1cs(18)/.00000000000000506E0/
  DATA atn1cs(19)/ - .00000000000000082E0/
  DATA atn1cs(20)/.00000000000000013E0/
  DATA atn1cs(21)/ - .00000000000000002E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  R9ATN1
  IF ( first ) THEN
    eps = R1MACH(3)
    ntatn1 = INITS(atn1cs,21,0.1*eps)
    !
    xsml = SQRT(0.1*eps)
    xbig = 1.571/SQRT(eps)
    xmax = 1.571/eps
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1.0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','R9ATN1',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',&
      2,2)
    IF ( y>xbig ) CALL XERMSG('SLATEC','R9ATN1',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'&
      ,1,1)
    !
    R9ATN1 = (ATAN(X)-X)/X**3
    GOTO 99999
  ENDIF
  !
  IF ( y<=xsml ) R9ATN1 = -1.0/3.0
  IF ( y<=xsml ) RETURN
  !
  R9ATN1 = -0.25 + CSEVL(2.0*y*y-1.,atn1cs,ntatn1)
  RETURN
  !
  99999 END FUNCTION R9ATN1
