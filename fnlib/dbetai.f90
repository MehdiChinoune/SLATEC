!DECK DBETAI
REAL(8) FUNCTION DBETAI(X,Pin,Qin)
  IMPLICIT NONE
  INTEGER i, ib, n
  !***BEGIN PROLOGUE  DBETAI
  !***PURPOSE  Calculate the incomplete Beta function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7F
  !***TYPE      DOUBLE PRECISION (BETAI-S, DBETAI-D)
  !***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !   DBETAI calculates the DOUBLE PRECISION incomplete beta function.
  !
  !   The incomplete beta function ratio is the probability that a
  !   random variable from a beta distribution having parameters PIN and
  !   QIN will be less than or equal to X.
  !
  !     -- Input Arguments -- All arguments are DOUBLE PRECISION.
  !   X      upper limit of integration.  X must be in (0,1) inclusive.
  !   PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
  !   QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
  !
  !***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
  !                 179, Communications of the ACM 17, 3 (March 1974),
  !                 pp. 156.
  !***ROUTINES CALLED  D1MACH, DLBETA, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  !***END PROLOGUE  DBETAI
  REAL(8) :: X, Pin, Qin, alneps, alnsml, c, eps, finsum, p, &
    ps, q, sml, term, xb, xi, y, D1MACH, DLBETA, p1
  LOGICAL first
  SAVE eps, alneps, sml, alnsml, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBETAI
  IF ( first ) THEN
    eps = D1MACH(3)
    alneps = LOG(eps)
    sml = D1MACH(1)
    alnsml = LOG(sml)
  ENDIF
  first = .FALSE.
  !
  IF ( X<0.D0.OR.X>1.D0 )&
    CALL XERMSG('SLATEC','DBETAI','X IS NOT IN THE RANGE (0,1)',1,2)
  IF ( Pin<=0.D0.OR.Qin<=0.D0 )&
    CALL XERMSG('SLATEC','DBETAI','P AND/OR Q IS LE ZERO',2,2)
  !
  y = X
  p = Pin
  q = Qin
  IF ( q>p.OR.X>=0.8D0 ) THEN
    IF ( X>=0.2D0 ) THEN
      y = 1.0D0 - y
      p = Qin
      q = Pin
    ENDIF
  ENDIF
  !
  IF ( (p+q)*y/(p+1.D0)<eps ) THEN
    !
    DBETAI = 0.0D0
    xb = p*LOG(MAX(y,sml)) - LOG(p) - DLBETA(p,q)
    IF ( xb>alnsml.AND.y/=0.0D0 ) DBETAI = EXP(xb)
    IF ( y/=X.OR.p/=Pin ) DBETAI = 1.0D0 - DBETAI
    RETURN
  ELSE
    !
    ! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
    ! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
    !
    ps = q - AINT(q)
    IF ( ps==0.D0 ) ps = 1.0D0
    xb = p*LOG(y) - DLBETA(ps,p) - LOG(p)
    DBETAI = 0.0D0
    IF ( xb>=alnsml ) THEN
      !
      DBETAI = EXP(xb)
      term = DBETAI*p
      IF ( ps/=1.0D0 ) THEN
        n = INT( MAX(alneps/LOG(y),4.0D0) )
        DO i = 1, n
          xi = i
          term = term*(xi-ps)*y/xi
          DBETAI = DBETAI + term/(p+xi)
        ENDDO
      ENDIF
    ENDIF
    !
    ! NOW EVALUATE THE FINITE SUM, MAYBE.
    !
    IF ( q>1.0D0 ) THEN
      !
      xb = p*LOG(y) + q*LOG(1.0D0-y) - DLBETA(p,q) - LOG(q)
      ib = INT( MAX(xb/alnsml,0.0D0) )
      term = EXP(xb-ib*alnsml)
      c = 1.0D0/(1.D0-y)
      p1 = q*c/(p+q-1.D0)
      !
      finsum = 0.0D0
      n = INT( q )
      IF ( q==REAL(n, 8) ) n = n - 1
      DO i = 1, n
        IF ( p1<=1.0D0.AND.term/eps<=finsum ) EXIT
        xi = i
        term = (q-xi+1.0D0)*c*term/(p+q-xi)
        !
        IF ( term>1.0D0 ) ib = ib - 1
        IF ( term>1.0D0 ) term = term*sml
        !
        IF ( ib==0 ) finsum = finsum + term
      ENDDO
      !
      DBETAI = DBETAI + finsum
    ENDIF
  ENDIF
  IF ( y/=X.OR.p/=Pin ) DBETAI = 1.0D0 - DBETAI
  DBETAI = MAX(MIN(DBETAI,1.0D0),0.0D0)
  RETURN
END FUNCTION DBETAI
