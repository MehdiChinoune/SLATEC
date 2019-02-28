!DECK BETAI
REAL FUNCTION BETAI(X,Pin,Qin)
  IMPLICIT NONE
  REAL ALBETA, alneps, alnsml, c, eps, finsum, p, p1, Pin, ps, q, &
    Qin, R1MACH, sml, term, X, xb, y
  INTEGER i, ib, n
  !***BEGIN PROLOGUE  BETAI
  !***PURPOSE  Calculate the incomplete Beta function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7F
  !***TYPE      SINGLE PRECISION (BETAI-S, DBETAI-D)
  !***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !   BETAI calculates the REAL incomplete beta function.
  !
  !   The incomplete beta function ratio is the probability that a
  !   random variable from a beta distribution having parameters PIN and
  !   QIN will be less than or equal to X.
  !
  !     -- Input Arguments -- All arguments are REAL.
  !   X      upper limit of integration.  X must be in (0,1) inclusive.
  !   PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
  !   QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
  !
  !***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
  !                 179, Communications of the ACM 17, 3 (March 1974),
  !                 pp. 156.
  !***ROUTINES CALLED  ALBETA, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  !***END PROLOGUE  BETAI
  LOGICAL first
  SAVE eps, alneps, sml, alnsml, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  BETAI
  IF ( first ) THEN
    eps = R1MACH(3)
    alneps = LOG(eps)
    sml = R1MACH(1)
    alnsml = LOG(sml)
  ENDIF
  first = .FALSE.
  !
  IF ( X<0..OR.X>1.0 ) CALL XERMSG('SLATEC','BETAI',&
    'X IS NOT IN THE RANGE (0,1)',1,2)
  IF ( Pin<=0..OR.Qin<=0. )&
    CALL XERMSG('SLATEC','BETAI','P AND/OR Q IS LE ZERO',2,2)
  !
  y = X
  p = Pin
  q = Qin
  IF ( q>p.OR.X>=0.8 ) THEN
    IF ( X>=0.2 ) THEN
      y = 1.0 - y
      p = Qin
      q = Pin
    ENDIF
  ENDIF
  !
  IF ( (p+q)*y/(p+1.)<eps ) THEN
    !
    BETAI = 0.0
    xb = p*LOG(MAX(y,sml)) - LOG(p) - ALBETA(p,q)
    IF ( xb>alnsml.AND.y/=0. ) BETAI = EXP(xb)
    IF ( y/=X.OR.p/=Pin ) BETAI = 1.0 - BETAI
    RETURN
  ELSE
    !
    ! EVALUATE THE INFINITE SUM FIRST.
    ! TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
    !
    ps = q - AINT(q)
    IF ( ps==0. ) ps = 1.0
    xb = p*LOG(y) - ALBETA(ps,p) - LOG(p)
    BETAI = 0.0
    IF ( xb>=alnsml ) THEN
      !
      BETAI = EXP(xb)
      term = BETAI*p
      IF ( ps/=1.0 ) THEN
        !
        n = INT( MAX(alneps/LOG(y),4.0E0) )
        DO i = 1, n
          term = term*(i-ps)*y/i
          BETAI = BETAI + term/(p+i)
        ENDDO
      ENDIF
    ENDIF
    !
    ! NOW EVALUATE THE FINITE SUM, MAYBE.
    !
    IF ( q>1.0 ) THEN
      !
      xb = p*LOG(y) + q*LOG(1.0-y) - ALBETA(p,q) - LOG(q)
      ib = INT( MAX(xb/alnsml,0.0E0) )
      term = EXP(xb-ib*alnsml)
      c = 1.0/(1.0-y)
      p1 = q*c/(p+q-1.)
      !
      finsum = 0.0
      n = INT( q )
      IF ( q==REAL(n) ) n = n - 1
      DO i = 1, n
        IF ( p1<=1.0.AND.term/eps<=finsum ) EXIT
        term = (q-i+1)*c*term/(p+q-i)
        !
        IF ( term>1.0 ) ib = ib - 1
        IF ( term>1.0 ) term = term*sml
        !
        IF ( ib==0 ) finsum = finsum + term
      ENDDO
      !
      BETAI = BETAI + finsum
    ENDIF
  ENDIF
  IF ( y/=X.OR.p/=Pin ) BETAI = 1.0 - BETAI
  BETAI = MAX(MIN(BETAI,1.0),0.0)
  RETURN
END FUNCTION BETAI
