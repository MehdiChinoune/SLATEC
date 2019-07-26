!** BETAI
REAL(SP) ELEMENTAL FUNCTION BETAI(X,Pin,Qin)
  !> Calculate the incomplete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7F
  !***
  ! **Type:**      SINGLE PRECISION (BETAI-S, DBETAI-D)
  !***
  ! **Keywords:**  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   BETAI calculates the REAL incomplete beta function.
  !
  !   The incomplete beta function ratio is the probability that a
  !   random variable from a beta distribution having parameters PIN and
  !   QIN will be less than or equal to X.
  !
  !     -- Input Arguments -- All arguments are REAL.
  !   X      upper limit of integration.  X must be in (0,1) inclusive.
  !   PIN    first beta distribution parameter.  PIN must be > 0.0.
  !   QIN    second beta distribution parameter.  QIN must be > 0.0.
  !
  !***
  ! **References:**  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
  !                 179, Communications of the ACM 17, 3 (March 1974),
  !                 pp. 156.
  !***
  ! **Routines called:**  ALBETA, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : eps_2_sp, tiny_sp
  !
  REAL(SP), INTENT(IN) :: X, Pin, Qin
  !
  INTEGER :: i, ib, n
  REAL(SP) :: c, finsum, p, p1, ps, q, term, xb, y
  REAL(SP), PARAMETER :: eps = eps_2_sp, alneps = LOG(eps), sml = tiny_sp, &
    alnsml = LOG(sml)
  !* FIRST EXECUTABLE STATEMENT  BETAI
  !
  IF( X<0. .OR. X>1._SP ) THEN
    ERROR STOP 'BETAI : X IS NOT IN THE RANGE (0,1)'
  ELSEIF( Pin<=0. .OR. Qin<=0. ) THEN
    ERROR STOP 'BETAI : P AND/OR Q IS <= 0'
  END IF
  !
  y = X
  p = Pin
  q = Qin
  IF( ( q>p .AND. X>=0.2 ) .OR. X>=0.8 ) THEN
    y = 1._SP - y
    p = Qin
    q = Pin
  END IF
  !
  IF( (p+q)*y/(p+1._SP)<eps ) THEN
    !
    BETAI = 0._SP
    xb = p*LOG(MAX(y,sml)) - LOG(p) - ALBETA(p,q)
    IF( xb>alnsml .AND. y/=0. ) BETAI = EXP(xb)
    IF( y/=X .OR. p/=Pin ) BETAI = 1._SP - BETAI
  ELSE
    !
    ! EVALUATE THE INFINITE SUM FIRST.
    ! TERM WILL EQUAL Y**P/BETA(PS,P) * (1.-PS)I * Y**I / FAC(I)
    !
    ps = q - AINT(q)
    IF( ps==0. ) ps = 1._SP
    xb = p*LOG(y) - ALBETA(ps,p) - LOG(p)
    BETAI = 0._SP
    IF( xb>=alnsml ) THEN
      !
      BETAI = EXP(xb)
      term = BETAI*p
      IF( ps/=1._SP ) THEN
        !
        n = INT( MAX(alneps/LOG(y),4._SP) )
        DO i = 1, n
          term = term*(i-ps)*y/i
          BETAI = BETAI + term/(p+i)
        END DO
      END IF
    END IF
    !
    ! NOW EVALUATE THE FINITE SUM, MAYBE.
    !
    IF( q>1._SP ) THEN
      !
      xb = p*LOG(y) + q*LOG(1._SP-y) - ALBETA(p,q) - LOG(q)
      ib = INT( MAX(xb/alnsml,0._SP) )
      term = EXP(xb-ib*alnsml)
      c = 1._SP/(1._SP-y)
      p1 = q*c/(p+q-1._SP)
      !
      finsum = 0._SP
      n = INT( q )
      IF( q==REAL(n) ) n = n - 1
      DO i = 1, n
        IF( p1<=1._SP .AND. term/eps<=finsum ) EXIT
        term = (q-i+1)*c*term/(p+q-i)
        !
        IF( term>1._SP ) ib = ib - 1
        IF( term>1._SP ) term = term*sml
        !
        IF( ib==0 ) finsum = finsum + term
      END DO
      !
      BETAI = BETAI + finsum
    END IF
    IF( y/=X .OR. p/=Pin ) BETAI = 1._SP - BETAI
    BETAI = MAX(MIN(BETAI,1._SP),0._SP)
  END IF

  RETURN
END FUNCTION BETAI