!** R9GMIC
REAL FUNCTION R9GMIC(A,X,Alx)
  !>
  !***
  !  Compute the complementary incomplete Gamma function for A
  !            near a negative integer and for small X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (R9GMIC-S, D9GMIC-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the complementary incomplete gamma function for A near
  ! a negative integer and for small X.
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

  REAL A, alng, Alx, fk, fkp1, fm, s, sgng, t, te, X
  INTEGER k, m, ma, mm1
  REAL, PARAMETER :: euler = .5772156649015329E0
  REAL :: eps = 0., bot = 0.
  !* FIRST EXECUTABLE STATEMENT  R9GMIC
  IF ( eps==0.0 ) eps = 0.5*R1MACH(3)
  IF ( bot==0.0 ) bot = LOG(R1MACH(1))
  !
  IF ( A>0.0 ) CALL XERMSG('SLATEC','R9GMIC',&
    'A MUST BE NEAR A NEGATIVE INTEGER',2,2)
  IF ( X<=0.0 ) CALL XERMSG('SLATEC','R9GMIC','X MUST BE GT ZERO',3,2)
  !
  ma = INT( A - 0.5 )
  fm = -ma
  m = -ma
  !
  te = 1.0
  t = 1.0
  s = t
  DO k = 1, 200
    fkp1 = k + 1
    te = -X*te/(fm+fkp1)
    t = te/fkp1
    s = s + t
    IF ( ABS(t)<eps*s ) GOTO 100
  END DO
  CALL XERMSG('SLATEC','R9GMIC',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',4,2)
  !
  100  R9GMIC = -Alx - euler + X*s/(fm+1.0)
  IF ( m==0 ) RETURN
  !
  IF ( m==1 ) R9GMIC = -R9GMIC - 1.0 + 1.0/X
  IF ( m==1 ) RETURN
  !
  te = fm
  t = 1.0
  s = t
  mm1 = m - 1
  DO k = 1, mm1
    fk = k
    te = -X*te/fk
    t = te/(fm-fk)
    s = s + t
    IF ( ABS(t)<eps*ABS(s) ) EXIT
  END DO
  !
  DO k = 1, m
    R9GMIC = R9GMIC + 1.0/k
  END DO
  !
  sgng = 1.0
  IF ( MOD(m,2)==1 ) sgng = -1.0
  alng = LOG(R9GMIC) - LOG_GAMMA(fm+1.0)
  !
  R9GMIC = 0.0
  IF ( alng>bot ) R9GMIC = sgng*EXP(alng)
  IF ( s/=0.0 ) R9GMIC = R9GMIC + SIGN(EXP(-fm*Alx+LOG(ABS(s)/fm)),s)
  !
  IF ( R9GMIC==0.0.AND.s==0.0 )&
    CALL XERMSG('SLATEC','R9GMIC','RESULT UNDERFLOWS',1,1)
  !
END FUNCTION R9GMIC
