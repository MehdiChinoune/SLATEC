!** GAMIC
REAL(SP) ELEMENTAL FUNCTION GAMIC(A,X)
  !> Calculate the complementary incomplete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (GAMIC-S, DGAMIC-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Evaluate the complementary incomplete gamma function
  !
  !   GAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
  !
  !   GAMIC is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though GAMIC is defined for X <
  !   0.0), except that for X = 0 and A <= 0.0, GAMIC is undefined.
  !
  !   GAMIC, A, and X are REAL.
  !
  !   A slight deterioration of 2 or 3 digits accuracy will occur when
  !   GAMIC is very large or very small in absolute value, because log-
  !   arithmic variables are used.  Also, if the parameter A is very close
  !   to a negative integer (but not a negative integer), there is a loss
  !   of accuracy, which is reported if the result is less than half
  !   machine precision.
  !
  !***
  ! **References:**  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***
  ! **Routines called:**  ALGAMS, R1MACH, R9GMIC, R9GMIT, R9LGIC,
  !                    R9LGIT, XERCLR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : eps_2_sp, eps_sp, tiny_sp
  !
  REAL(SP), INTENT(IN) :: A, X
  !
  INTEGER :: izero, ma
  REAL(SP) :: aeps, algap1, alngs, alx, e, fm, gstar, h, sga, sgng,  sgngam, sgngs, t
  REAL(SP), PARAMETER :: eps = 0.5_SP*eps_2_sp, sqeps = SQRT(eps_sp), &
    alneps = -LOG(eps_2_sp), bot = LOG(tiny_sp)
  !* FIRST EXECUTABLE STATEMENT  GAMIC
  !
  IF( X<0._SP ) THEN
    ERROR STOP 'GAMIC : X IS NEGATIVE'
  ELSEIF( X>0._SP ) THEN
    !
    alx = LOG(X)
    sga = 1._SP
    IF( A/=0._SP ) sga = SIGN(1._SP,A)
    ma = INT( A + 0.5_SP*sga )
    aeps = A - ma
    !
    izero = 0
    IF( X>=1._SP ) THEN
      !
      IF( A<X ) THEN
        GAMIC = EXP(R9LGIC(A,X,alx))
        RETURN
      END IF
      !
      sgngam = 1._SP
      algap1 = LOG_GAMMA(A+1._SP)
      sgngs = 1._SP
      alngs = R9LGIT(A,X,algap1)
    ELSE
      !
      IF( A<=0.5_SP .AND. ABS(aeps)<=0.001 ) THEN
        fm = -ma
        e = 2._SP
        IF( fm>1._SP ) e = 2._SP*(fm+2._SP)/(fm*fm-1._SP)
        e = e - alx*X**(-0.001_SP)
        IF( e*ABS(aeps)<=eps ) THEN
          !
          GAMIC = R9GMIC(A,X,alx)
          RETURN
        END IF
      END IF
      !
      CALL ALGAMS(A+1._SP,algap1,sgngam)
      gstar = R9GMIT(A,X,algap1,sgngam)
      IF( gstar==0._SP ) izero = 1
      IF( gstar/=0._SP ) alngs = LOG(ABS(gstar))
      IF( gstar/=0._SP ) sgngs = SIGN(1._SP,gstar)
    END IF
    !
    ! EVALUATION OF GAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
    !
    h = 1._SP
    IF( izero/=1 ) THEN
      !
      t = A*alx + alngs
      IF( t>alneps ) THEN
        !
        sgng = -sgngs*sga*sgngam
        t = t + algap1 - LOG(ABS(A))
        GAMIC = sgng*EXP(t)
        RETURN
      ELSEIF( t>(-alneps) ) THEN
        h = 1._SP - sgngs*EXP(t)
        ! IF( ABS(h)<sqeps ) CALL XERMSG('GAMIC','RESULT LT HALF PRECISION',1,1)
      END IF
    END IF
  ELSE
    IF( A<=0._SP ) ERROR STOP 'GAMIC : X = 0 AND A LE 0 SO GAMIC IS UNDEFINED'
    !
    GAMIC = EXP(LOG_GAMMA(A+1._SP)-LOG(A))
    RETURN
  END IF
  !
  sgng = SIGN(1._SP,h)*sga*sgngam
  t = LOG(ABS(h)) + algap1 - LOG(ABS(A))
  GAMIC = sgng*EXP(t)

  RETURN
END FUNCTION GAMIC