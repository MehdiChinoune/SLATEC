!*==DBINOM.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBINOM
DOUBLE PRECISION FUNCTION DBINOM(N,M)
  IMPLICIT NONE
  !*--DBINOM5
  !*** Start of declarations inserted by SPAG
  INTEGER i , k , M , N
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBINOM
  !***PURPOSE  Compute the binomial coefficients.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C1
  !***TYPE      DOUBLE PRECISION (BINOM-S, DBINOM-D)
  !***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBINOM(N,M) calculates the double precision binomial coefficient
  ! for integer arguments N and M.  The result is (N!)/((M!)(N-M)!).
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9LGMC, DLNREL, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBINOM
  DOUBLE PRECISION corr , fintmx , sq2pil , xk , xn , xnk , D9LGMC , &
    DLNREL , D1MACH , bilnmx
  LOGICAL first
  SAVE sq2pil , bilnmx , fintmx , first
  DATA sq2pil/0.91893853320467274178032973640562D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBINOM
  IF ( first ) THEN
    bilnmx = LOG(D1MACH(2)) - 0.0001D0
    fintmx = 0.9D0/D1MACH(3)
  ENDIF
  first = .FALSE.
  !
  IF ( N<0.OR.M<0 ) CALL XERMSG('SLATEC','DBINOM','N OR M LT ZERO',1,2)
  IF ( N<M ) CALL XERMSG('SLATEC','DBINOM','N LT M',2,2)
  !
  k = MIN(M,N-M)
  IF ( k<=20 ) THEN
    IF ( k*LOG(AMAX0(N,1))<=bilnmx ) THEN
      !
      DBINOM = 1.0D0
      IF ( k==0 ) RETURN
      DO i = 1 , k
        xn = N - i + 1
        xk = i
        DBINOM = DBINOM*(xn/xk)
      ENDDO
      !
      IF ( DBINOM<fintmx ) DBINOM = AINT(DBINOM+0.5D0)
      RETURN
    ENDIF
  ENDIF
  !
  ! IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
  IF ( k<9 ) CALL XERMSG('SLATEC','DBINOM',&
    'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG',3,2)
  !
  xn = N + 1
  xk = k + 1
  xnk = N - k + 1
  !
  corr = D9LGMC(xn) - D9LGMC(xk) - D9LGMC(xnk)
  DBINOM = xk*LOG(xnk/xk) - xn*DLNREL(-(xk-1.0D0)/xn) - 0.5D0*LOG(xn*xnk/xk)&
    + 1.0D0 - sq2pil + corr
  !
  IF ( DBINOM>bilnmx ) CALL XERMSG('SLATEC','DBINOM',&
    'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG'&
    ,3,2)
  !
  DBINOM = EXP(DBINOM)
  IF ( DBINOM<fintmx ) DBINOM = AINT(DBINOM+0.5D0)
  !
END FUNCTION DBINOM
