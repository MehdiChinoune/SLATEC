!** D9GMIC
REAL(DP) FUNCTION D9GMIC(A,X,Alx)
  !> Compute the complementary incomplete Gamma function for A
  !            near a negative integer and X small.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (R9GMIC-S, D9GMIC-D)
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
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  REAL(DP) :: A, X, Alx
  INTEGER :: k, m, mm1
  REAL(DP) :: alng, fk, fkp1, fm, s, sgng, t, te
  REAL(DP), PARAMETER :: eps = 0.5_DP*D1MACH(3), bot = LOG(D1MACH(1))
  REAL(DP), PARAMETER :: euler = 0.57721566490153286060651209008240_DP
  !* FIRST EXECUTABLE STATEMENT  D9GMIC
  !
  IF( A>0._DP ) CALL XERMSG('D9GMIC','A MUST BE NEAR A NEGATIVE INTEGER',2,2)
  IF( X<=0._DP ) CALL XERMSG('D9GMIC','X MUST BE GT ZERO',3,2)
  !
  m = INT( -(A-0.5_DP) )
  fm = m
  !
  te = 1._DP
  t = 1._DP
  s = t
  DO k = 1, 200
    fkp1 = k + 1
    te = -X*te/(fm+fkp1)
    t = te/fkp1
    s = s + t
    IF( ABS(t)<eps*s ) GOTO 100
  END DO
  CALL XERMSG('D9GMIC',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',4,2)
  !
  100  D9GMIC = -Alx - euler + X*s/(fm+1._DP)
  IF( m==0 ) RETURN
  !
  IF( m==1 ) D9GMIC = -D9GMIC - 1._DP + 1._DP/X
  IF( m==1 ) RETURN
  !
  te = fm
  t = 1._DP
  s = t
  mm1 = m - 1
  DO k = 1, mm1
    fk = k
    te = -X*te/fk
    t = te/(fm-fk)
    s = s + t
    IF( ABS(t)<eps*ABS(s) ) EXIT
  END DO
  !
  DO k = 1, m
    D9GMIC = D9GMIC + 1._DP/k
  END DO
  !
  sgng = 1._DP
  IF( MOD(m,2)==1 ) sgng = -1._DP
  alng = LOG(D9GMIC) - LOG_GAMMA(fm+1._DP)
  !
  D9GMIC = 0._DP
  IF( alng>bot ) D9GMIC = sgng*EXP(alng)
  IF( s/=0._DP ) D9GMIC = D9GMIC + SIGN(EXP(-fm*Alx+LOG(ABS(s)/fm)),s)
  !
  IF( D9GMIC==0._DP .AND. s==0._DP )&
    CALL XERMSG('D9GMIC','RESULT UNDERFLOWS',1,1)
  !
END FUNCTION D9GMIC
