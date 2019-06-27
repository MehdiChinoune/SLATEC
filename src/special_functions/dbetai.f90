!** DBETAI
REAL(DP) ELEMENTAL FUNCTION DBETAI(X,Pin,Qin)
  !> Calculate the incomplete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7F
  !***
  ! **Type:**      DOUBLE PRECISION (BETAI-S, DBETAI-D)
  !***
  ! **Keywords:**  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   DBETAI calculates the DOUBLE PRECISION incomplete beta function.
  !
  !   The incomplete beta function ratio is the probability that a
  !   random variable from a beta distribution having parameters PIN and
  !   QIN will be less than or equal to X.
  !
  !     -- Input Arguments -- All arguments are DOUBLE PRECISION.
  !   X      upper limit of integration.  X must be in (0,1) inclusive.
  !   PIN    first beta distribution parameter.  PIN must be > 0.0.
  !   QIN    second beta distribution parameter.  QIN must be > 0.0.
  !
  !***
  ! **References:**  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
  !                 179, Communications of the ACM 17, 3 (March 1974), pp. 156.
  !***
  ! **Routines called:**  D1MACH, DLBETA, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: X, Pin, Qin
  INTEGER :: i, ib, n
  REAL(DP) :: c, finsum, p, ps, q, term, xb, xi, y, p1
  REAL(DP), PARAMETER :: eps = D1MACH(3), alneps = LOG(eps), sml = D1MACH(1), &
    alnsml = LOG(sml)
  !* FIRST EXECUTABLE STATEMENT  DBETAI
  !
  IF( X<0._DP .OR. X>1._DP ) THEN
    ERROR STOP 'DBETAI : X IS NOT IN THE RANGE (0,1)'
  ELSEIF( Pin<=0._DP .OR. Qin<=0._DP ) THEN
    ERROR STOP 'DBETAI : P AND/OR Q IS <= 0'
  END IF
  !
  y = X
  p = Pin
  q = Qin
  IF( q>p .OR. X>=0.8_DP ) THEN
    IF( X>=0.2_DP ) THEN
      y = 1._DP - y
      p = Qin
      q = Pin
    END IF
  END IF
  !
  IF( (p+q)*y/(p+1._DP)<eps ) THEN
    !
    DBETAI = 0._DP
    xb = p*LOG(MAX(y,sml)) - LOG(p) - DLBETA(p,q)
    IF( xb>alnsml .AND. y/=0._DP ) DBETAI = EXP(xb)
    IF( y/=X .OR. p/=Pin ) DBETAI = 1._DP - DBETAI
  ELSE
    !
    ! EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
    ! Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
    !
    ps = q - AINT(q)
    IF( ps==0._DP ) ps = 1._DP
    xb = p*LOG(y) - DLBETA(ps,p) - LOG(p)
    DBETAI = 0._DP
    IF( xb>=alnsml ) THEN
      !
      DBETAI = EXP(xb)
      term = DBETAI*p
      IF( ps/=1._DP ) THEN
        n = INT( MAX(alneps/LOG(y),4._DP) )
        DO i = 1, n
          xi = i
          term = term*(xi-ps)*y/xi
          DBETAI = DBETAI + term/(p+xi)
        END DO
      END IF
    END IF
    !
    ! NOW EVALUATE THE FINITE SUM, MAYBE.
    !
    IF( q>1._DP ) THEN
      !
      xb = p*LOG(y) + q*LOG(1._DP-y) - DLBETA(p,q) - LOG(q)
      ib = INT( MAX(xb/alnsml,0._DP) )
      term = EXP(xb-ib*alnsml)
      c = 1._DP/(1._DP-y)
      p1 = q*c/(p+q-1._DP)
      !
      finsum = 0._DP
      n = INT( q )
      IF( q==REAL( n, DP ) ) n = n - 1
      DO i = 1, n
        IF( p1<=1._DP .AND. term/eps<=finsum ) EXIT
        xi = i
        term = (q-xi+1._DP)*c*term/(p+q-xi)
        !
        IF( term>1._DP ) ib = ib - 1
        IF( term>1._DP ) term = term*sml
        !
        IF( ib==0 ) finsum = finsum + term
      END DO
      !
      DBETAI = DBETAI + finsum
    END IF
    IF( y/=X .OR. p/=Pin ) DBETAI = 1._DP - DBETAI
    DBETAI = MAX(MIN(DBETAI,1._DP),0._DP)
  END IF

  RETURN
END FUNCTION DBETAI