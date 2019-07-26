!** DGAMIC
REAL(DP) ELEMENTAL FUNCTION DGAMIC(A,X)
  !> Calculate the complementary incomplete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (GAMIC-S, DGAMIC-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Evaluate the complementary incomplete Gamma function
  !
  !   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
  !
  !   DGAMIC is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though DGAMIC is defined for X <
  !   0.0), except that for X = 0 and A <= 0.0, DGAMIC is undefined.
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
  !***
  ! **References:**  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***
  ! **Routines called:**  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS,
  !                    DLNGAM, XERCLR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : eps_2_dp, eps_dp, tiny_dp
  !
  REAL(DP), INTENT(IN) :: A, X
  !
  INTEGER :: izero
  REAL(DP) :: aeps, ainta, algap1, alngs, alx, e, gstar, h, sga, sgng, sgngam, sgngs, t
  REAL(DP), PARAMETER :: eps = 0.5_DP*eps_2_dp, sqeps = SQRT(eps_dp), &
    alneps = -LOG(eps_2_dp), bot = LOG(tiny_dp)
  !* FIRST EXECUTABLE STATEMENT  DGAMIC
  !
  IF( X<0._DP ) ERROR STOP 'DGAMIC : X IS NEGATIVE'
  !
  IF( X>0._DP ) THEN
    !
    alx = LOG(X)
    sga = 1._DP
    IF( A/=0._DP ) sga = SIGN(1._DP,A)
    ainta = AINT(A+0.5_DP*sga)
    aeps = A - ainta
    !
    izero = 0
    IF( X>=1._DP ) THEN
      !
      IF( A<X ) THEN
        DGAMIC = EXP(D9LGIC(A,X,alx))
        RETURN
      END IF
      !
      sgngam = 1._DP
      algap1 = LOG_GAMMA(A+1._DP)
      sgngs = 1._DP
      alngs = D9LGIT(A,X,algap1)
    ELSE
      !
      IF( A<=0.5_DP .AND. ABS(aeps)<=0.001_DP ) THEN
        e = 2._DP
        IF( -ainta>1._DP ) e = 2._DP*(-ainta+2._DP)/(ainta*ainta-1._DP)
        e = e - alx*X**(-0.001_DP)
        IF( e*ABS(aeps)<=eps ) THEN
          !
          DGAMIC = D9GMIC(A,X,alx)
          RETURN
        END IF
      END IF
      !
      CALL DLGAMS(A+1._DP,algap1,sgngam)
      gstar = D9GMIT(A,X,algap1,sgngam)
      IF( gstar==0._DP ) izero = 1
      IF( gstar/=0._DP ) alngs = LOG(ABS(gstar))
      IF( gstar/=0._DP ) sgngs = SIGN(1._DP,gstar)
    END IF
    !
    ! EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
    !
    h = 1._DP
    IF( izero/=1 ) THEN
      t = A*alx + alngs
      IF( t>alneps ) THEN
        sgng = -sgngs*sga*sgngam
        t = t + algap1 - LOG(ABS(A))
        DGAMIC = sgng*EXP(t)
        RETURN
      ELSE
        IF( t>(-alneps) ) h = 1._DP - sgngs*EXP(t)
        ! IF( ABS(h)<sqeps ) CALL XERMSG('DGAMIC','RESULT LT HALF PRECISION',1,1)
      END IF
    END IF
    sgng = SIGN(1._DP,h)*sga*sgngam
    t = LOG(ABS(h)) + algap1 - LOG(ABS(A))
    DGAMIC = sgng*EXP(t)
  ELSEIF( A<=0._DP ) THEN
    ERROR STOP 'DGAMIC : X = 0 AND A <= 0 SO DGAMIC IS UNDEFINED'
  ELSE
    DGAMIC = EXP(LOG_GAMMA(A+1._DP)-LOG(A))
  END IF

  RETURN
END FUNCTION DGAMIC