!DECK DGAMIC
REAL(8) FUNCTION DGAMIC(A,X)
  IMPLICIT NONE
  INTEGER izero
  !***BEGIN PROLOGUE  DGAMIC
  !***PURPOSE  Calculate the complementary incomplete Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      DOUBLE PRECISION (GAMIC-S, DGAMIC-D)
  !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !   Evaluate the complementary incomplete Gamma function
  !
  !   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
  !
  !   DGAMIC is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though DGAMIC is defined for X .LT.
  !   0.0), except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.
  !
  !   DGAMIC, A, and X are DOUBLE PRECISION.
  !
  !   A slight deterioration of 2 or 3 digits accuracy will occur when
  !   DGAMIC is very large or very small in absolute value, because log-
  !   arithmic variables are used.  Also, if the parameter A is very close
  !   to a negative INTEGER (but not a negative integer), there is a loss
  !   of accuracy, which is reported if the result is less than half
  !   machine precision.
  !
  !***REFERENCES  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***ROUTINES CALLED  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS,
  !                    DLNGAM, XERCLR, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  !***END PROLOGUE  DGAMIC
  REAL(8) :: A, X, aeps, ainta, algap1, alneps, alngs, alx, &
    bot, e, eps, gstar, h, sga, sgng, sgngam, sgngs, &
    sqeps, t, D1MACH, DLNGAM, D9GMIC, D9GMIT, D9LGIC, &
    D9LGIT
  LOGICAL first
  SAVE eps, sqeps, alneps, bot, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DGAMIC
  IF ( first ) THEN
    eps = 0.5D0*D1MACH(3)
    sqeps = SQRT(D1MACH(4))
    alneps = -LOG(D1MACH(3))
    bot = LOG(D1MACH(1))
  ENDIF
  first = .FALSE.
  !
  IF ( X<0.D0 ) CALL XERMSG('SLATEC','DGAMIC','X IS NEGATIVE',2,2)
  !
  IF ( X>0.D0 ) THEN
    !
    alx = LOG(X)
    sga = 1.0D0
    IF ( A/=0.D0 ) sga = SIGN(1.0D0,A)
    ainta = AINT(A+0.5D0*sga)
    aeps = A - ainta
    !
    izero = 0
    IF ( X>=1.0D0 ) THEN
      !
      IF ( A<X ) THEN
        DGAMIC = EXP(D9LGIC(A,X,alx))
        RETURN
      ENDIF
      !
      sgngam = 1.0D0
      algap1 = DLNGAM(A+1.0D0)
      sgngs = 1.0D0
      alngs = D9LGIT(A,X,algap1)
    ELSE
      !
      IF ( A<=0.5D0.AND.ABS(aeps)<=0.001D0 ) THEN
        e = 2.0D0
        IF ( -ainta>1.D0 ) e = 2.D0*(-ainta+2.D0)/(ainta*ainta-1.0D0)
        e = e - alx*X**(-0.001D0)
        IF ( e*ABS(aeps)<=eps ) THEN
          !
          DGAMIC = D9GMIC(A,X,alx)
          RETURN
        ENDIF
      ENDIF
      !
      CALL DLGAMS(A+1.0D0,algap1,sgngam)
      gstar = D9GMIT(A,X,algap1,sgngam,alx)
      IF ( gstar==0.D0 ) izero = 1
      IF ( gstar/=0.D0 ) alngs = LOG(ABS(gstar))
      IF ( gstar/=0.D0 ) sgngs = SIGN(1.0D0,gstar)
    ENDIF
    !
    ! EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
    !
    h = 1.D0
    IF ( izero/=1 ) THEN
      !
      t = A*alx + alngs
      IF ( t>alneps ) THEN
        !
        sgng = -sgngs*sga*sgngam
        t = t + algap1 - LOG(ABS(A))
        IF ( t<bot ) CALL XERCLR
        DGAMIC = sgng*EXP(t)
        RETURN
      ELSE
        IF ( t>(-alneps) ) h = 1.0D0 - sgngs*EXP(t)
        !
        IF ( ABS(h)<sqeps ) CALL XERCLR
        IF ( ABS(h)<sqeps )&
          CALL XERMSG('SLATEC','DGAMIC','RESULT LT HALF PRECISION',1,1)
      ENDIF
    ENDIF
  ELSE
    IF ( A<=0.D0 ) CALL XERMSG('SLATEC','DGAMIC',&
      'X = 0 AND A LE 0 SO DGAMIC IS UNDEFINED',3,&
      2)
    !
    DGAMIC = EXP(DLNGAM(A+1.D0)-LOG(A))
    RETURN
  ENDIF
  !
  sgng = SIGN(1.0D0,h)*sga*sgngam
  t = LOG(ABS(h)) + algap1 - LOG(ABS(A))
  IF ( t<bot ) CALL XERCLR
  DGAMIC = sgng*EXP(t)
  RETURN
END FUNCTION DGAMIC
