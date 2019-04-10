!** DCHU
REAL(8) FUNCTION DCHU(A,B,X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithmic confluent hypergeometric function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C11
  !***
  ! **Type:**      DOUBLE PRECISION (CHU-S, DCHU-D)
  !***
  ! **Keywords:**  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DCHU(A,B,X) calculates the double precision logarithmic confluent
  ! hypergeometric function U(A,B,X) for double precision arguments
  ! A, B, and X.
  !
  ! This routine is not valid when 1+A-B is close to zero if X is small.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9CHU, DEXPRL, DGAMMA, DGAMR, DPOCH,
  !                    DPOCH1, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)

  INTEGER i, istrt, m, n
  REAL(8) :: A, B, X, aintb, alnx, a0, beps, b0, c0, factor, gamri1, gamrni, &
    pch1ai, pch1i, pochai, summ, t, xeps1, xi, xi1, xn, xtoeps
  REAL(8), EXTERNAL :: DGAMMA, D1MACH, DPOCH, DGAMR, DPOCH1, DEXPRL, D9CHU
  REAL(8), PARAMETER :: pi = 3.141592653589793238462643383279503D0
  REAL(8) :: eps = 0.0D0
  !* FIRST EXECUTABLE STATEMENT  DCHU
  IF ( eps==0.0D0 ) eps = D1MACH(3)
  !
  IF ( X==0.0D0 ) CALL XERMSG('SLATEC','DCHU',&
    'X IS ZERO SO DCHU IS INFINITE',1,2)
  IF ( X<0.0D0 ) CALL XERMSG('SLATEC','DCHU','X IS NEGATIVE, USE CCHU',2,2)
  !
  IF ( MAX(ABS(A),1.0D0)*MAX(ABS(1.0D0+A-B),1.0D0)>=0.99D0*ABS(X) ) THEN
    !
    ! THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
    ! APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
    !
    IF ( ABS(1.0D0+A-B)<SQRT(eps) ) CALL XERMSG('SLATEC','DCHU',&
      'ALGORITHMIS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X',10,2)
    !
    IF ( B>=0.0D0 ) aintb = AINT(B+0.5D0)
    IF ( B<0.0D0 ) aintb = AINT(B-0.5D0)
    beps = B - aintb
    n = INT( aintb )
    !
    alnx = LOG(X)
    xtoeps = EXP(-beps*alnx)
    !
    ! EVALUATE THE FINITE SUM.     -----------------------------------------
    !
    IF ( n>=1 ) THEN
      !
      ! NOW CONSIDER THE CASE B .GE. 1.0.
      !
      summ = 0.0D0
      m = n - 2
      IF ( m>=0 ) THEN
        t = 1.0D0
        summ = 1.0D0
        IF ( m/=0 ) THEN
          !
          DO i = 1, m
            xi = i
            t = t*(A-B+xi)*X/((1.0D0-B+xi)*xi)
            summ = summ + t
          END DO
        END IF
        !
        summ = DGAMMA(B-1.0D0)*DGAMR(A)*X**(1-n)*xtoeps*summ
      END IF
    ELSE
      !
      ! CONSIDER THE CASE B .LT. 1.0 FIRST.
      !
      summ = 1.0D0
      IF ( n/=0 ) THEN
        !
        t = 1.0D0
        m = -n
        DO i = 1, m
          xi1 = i - 1
          t = t*(A+xi1)*X/((B+xi1)*(xi1+1.0D0))
          summ = summ + t
        END DO
      END IF
      !
      summ = DPOCH(1.0D0+A-B,-A)*summ
    END IF
    !
    ! NEXT EVALUATE THE INFINITE SUM.     ----------------------------------
    !
    istrt = 0
    IF ( n<1 ) istrt = 1 - n
    xi = istrt
    !
    factor = (-1.0D0)**n*DGAMR(1.0D0+A-B)*X**istrt
    IF ( beps/=0.0D0 ) factor = factor*beps*pi/SIN(beps*pi)
    !
    pochai = DPOCH(A,xi)
    gamri1 = DGAMR(xi+1.0D0)
    gamrni = DGAMR(aintb+xi)
    b0 = factor*DPOCH(A,xi-beps)*gamrni*DGAMR(xi+1.0D0-beps)
    !
    IF ( ABS(xtoeps-1.0D0)<=0.5D0 ) THEN
      !
      ! X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE CAREFUL IN EVALUATING THE
      ! DIFFERENCES.
      !
      pch1ai = DPOCH1(A+xi,-beps)
      pch1i = DPOCH1(xi+1.0D0-beps,beps)
      c0 = factor*pochai*gamrni*gamri1*(-DPOCH1(B+xi,-beps)+pch1ai-pch1i+&
        beps*pch1ai*pch1i)
      !
      ! XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS)
      xeps1 = alnx*DEXPRL(-beps*alnx)
      !
      DCHU = summ + c0 + xeps1*b0
      xn = n
      DO i = 1, 1000
        xi = istrt + i
        xi1 = istrt + i - 1
        b0 = (A+xi1-beps)*b0*X/((xn+xi1)*(xi-beps))
        c0 = (A+xi1)*c0*X/((B+xi1)*xi)&
          - ((A-1.0D0)*(xn+2.D0*xi-1.0D0)+xi*(xi-beps))&
          *b0/(xi*(B+xi1)*(A+xi1-beps))
        t = c0 + xeps1*b0
        DCHU = DCHU + t
        IF ( ABS(t)<eps*ABS(DCHU) ) RETURN
      END DO
      CALL XERMSG('SLATEC','DCHU',&
        'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES',3,2)
    END IF
    !
    ! X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
    ! FORMULATION IS STABLE.
    !
    a0 = factor*pochai*DGAMR(B+xi)*gamri1/beps
    b0 = xtoeps*b0/beps
    !
    DCHU = summ + a0 - b0
    DO i = 1, 1000
      xi = istrt + i
      xi1 = istrt + i - 1
      a0 = (A+xi1)*a0*X/((B+xi1)*xi)
      b0 = (A+xi1-beps)*b0*X/((aintb+xi1)*(xi-beps))
      t = a0 - b0
      DCHU = DCHU + t
      IF ( ABS(t)<eps*ABS(DCHU) ) RETURN
    END DO
    CALL XERMSG('SLATEC','DCHU',&
      'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES',3,2)
  END IF
  !
  ! USE LUKE-S RATIONAL APPROXIMATION IN THE ASYMPTOTIC REGION.
  !
  DCHU = X**(-A)*D9CHU(A,B,X)
  !
  RETURN
END FUNCTION DCHU
