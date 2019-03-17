!DECK D9GMIC
REAL(8) FUNCTION D9GMIC(A,X,Alx)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  D9GMIC
  !***SUBSIDIARY
  !***PURPOSE  Compute the complementary incomplete Gamma function for A
  !            near a negative integer and X small.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      DOUBLE PRECISION (R9GMIC-S, D9GMIC-D)
  !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the complementary incomplete gamma function for A near
  ! a negative integer and for small X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  D9GMIC
  INTEGER k, m, mm1
  REAL(8) :: A, X, Alx, alng, bot, eps, euler, fk, fkp1, fm, &
    s, sgng, t, te, D1MACH, DLNGAM
  LOGICAL first
  SAVE euler, eps, bot, first
  DATA euler/0.57721566490153286060651209008240D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9GMIC
  IF ( first ) THEN
    eps = 0.5D0*D1MACH(3)
    bot = LOG(D1MACH(1))
  ENDIF
  first = .FALSE.
  !
  IF ( A>0.D0 ) CALL XERMSG('SLATEC','D9GMIC',&
    'A MUST BE NEAR A NEGATIVE INTEGER',2,2)
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','D9GMIC','X MUST BE GT ZERO',3,2)
  !
  m = INT( -(A-0.5D0) )
  fm = m
  !
  te = 1.0D0
  t = 1.0D0
  s = t
  DO k = 1, 200
    fkp1 = k + 1
    te = -X*te/(fm+fkp1)
    t = te/fkp1
    s = s + t
    IF ( ABS(t)<eps*s ) GOTO 100
  ENDDO
  CALL XERMSG('SLATEC','D9GMIC',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',4,2)
  !
  100  D9GMIC = -Alx - euler + X*s/(fm+1.0D0)
  IF ( m==0 ) RETURN
  !
  IF ( m==1 ) D9GMIC = -D9GMIC - 1.D0 + 1.D0/X
  IF ( m==1 ) RETURN
  !
  te = fm
  t = 1.D0
  s = t
  mm1 = m - 1
  DO k = 1, mm1
    fk = k
    te = -X*te/fk
    t = te/(fm-fk)
    s = s + t
    IF ( ABS(t)<eps*ABS(s) ) EXIT
  ENDDO
  !
  DO k = 1, m
    D9GMIC = D9GMIC + 1.0D0/k
  ENDDO
  !
  sgng = 1.0D0
  IF ( MOD(m,2)==1 ) sgng = -1.0D0
  alng = LOG(D9GMIC) - DLNGAM(fm+1.D0)
  !
  D9GMIC = 0.D0
  IF ( alng>bot ) D9GMIC = sgng*EXP(alng)
  IF ( s/=0.D0 ) D9GMIC = D9GMIC + SIGN(EXP(-fm*Alx+LOG(ABS(s)/fm)),s)
  !
  IF ( D9GMIC==0.D0.AND.s==0.D0 )&
    CALL XERMSG('SLATEC','D9GMIC','RESULT UNDERFLOWS',1,1)
  !
END FUNCTION D9GMIC
