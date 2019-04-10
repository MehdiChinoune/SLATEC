!** CHU
REAL FUNCTION CHU(A,B,X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithmic confluent hypergeometric function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C11
  !***
  ! **Type:**      SINGLE PRECISION (CHU-S, DCHU-D)
  !***
  ! **Keywords:**  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CHU computes the logarithmic confluent hypergeometric function,
  ! U(A,B,X).
  !
  ! Input Parameters:
  !       A   real
  !       B   real
  !       X   real and positive
  !
  ! This routine is not valid when 1+A-B is close to zero if X is small.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  EXPREL, GAMMA, GAMR, POCH, POCH1, R1MACH, R9CHU,
  !                    XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)

  REAL A, a0, aintb, alnx, B, b0, beps, c0, factor,gamri1, gamrni, pch1ai, &
    pch1i, pochai, summ, t, X, xeps1, xi, xi1, xn, xtoeps
  INTEGER i, istrt, m, n
  REAL, EXTERNAL :: EXPREL, GAMR, POCH1, R1MACH, R9CHU, GAMMA, POCH
  REAL, PARAMETER :: pi = 3.14159265358979324E0
  REAL :: eps = 0.0
  !* FIRST EXECUTABLE STATEMENT  CHU
  IF ( eps==0.0 ) eps = R1MACH(3)
  !
  IF ( X==0.0 ) CALL XERMSG('SLATEC','CHU','X IS ZERO SO CHU IS INFINITE',1,2)
  IF ( X<0.0 ) CALL XERMSG('SLATEC','CHU','X IS NEGATIVE, USE CCHU',2,2)
  !
  IF ( MAX(ABS(A),1.0)*MAX(ABS(1.0+A-B),1.0)>=0.99*ABS(X) ) THEN
    !
    ! THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
    ! APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
    !
    IF ( ABS(1.0+A-B)<SQRT(eps) ) CALL XERMSG('SLATEC','CHU',&
      'ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X',10,2)
    !
    aintb = AINT(B+0.5)
    IF ( B<0.0 ) aintb = AINT(B-0.5)
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
      summ = 0.0
      m = n - 2
      IF ( m>=0 ) THEN
        t = 1.0
        summ = 1.0
        IF ( m/=0 ) THEN
          !
          DO i = 1, m
            xi = i
            t = t*(A-B+xi)*X/((1.0-B+xi)*xi)
            summ = summ + t
          END DO
        END IF
        !
        summ = GAMMA(B-1.0)*GAMR(A)*X**(1-n)*xtoeps*summ
      END IF
    ELSE
      !
      ! CONSIDER THE CASE B .LT. 1.0 FIRST.
      !
      summ = 1.0
      IF ( n/=0 ) THEN
        !
        t = 1.0
        m = -n
        DO i = 1, m
          xi1 = i - 1
          t = t*(A+xi1)*X/((B+xi1)*(xi1+1.0))
          summ = summ + t
        END DO
      END IF
      !
      summ = POCH(1.0+A-B,-A)*summ
    END IF
    !
    ! NOW EVALUATE THE INFINITE SUM.     -----------------------------------
    !
    istrt = 0
    IF ( n<1 ) istrt = 1 - n
    xi = istrt
    !
    factor = (-1.0)**n*GAMR(1.0+A-B)*X**istrt
    IF ( beps/=0.0 ) factor = factor*beps*pi/SIN(beps*pi)
    !
    pochai = POCH(A,xi)
    gamri1 = GAMR(xi+1.0)
    gamrni = GAMR(aintb+xi)
    b0 = factor*POCH(A,xi-beps)*gamrni*GAMR(xi+1.0-beps)
    !
    IF ( ABS(xtoeps-1.0)<=0.5 ) THEN
      !
      ! X**(-BEPS) IS CLOSE TO 1.0, SO WE MUST BE CAREFUL IN EVALUATING
      ! THE DIFFERENCES
      !
      pch1ai = POCH1(A+xi,-beps)
      pch1i = POCH1(xi+1.0-beps,beps)
      c0 = factor*pochai*gamrni*gamri1*(-POCH1(B+xi,-beps)+pch1ai-pch1i+&
        beps*pch1ai*pch1i)
      !
      ! XEPS1 = (1.0 - X**(-BEPS)) / BEPS
      xeps1 = alnx*EXPREL(-beps*alnx)
      !
      CHU = summ + c0 + xeps1*b0
      xn = n
      DO i = 1, 1000
        xi = istrt + i
        xi1 = istrt + i - 1
        b0 = (A+xi1-beps)*b0*X/((xn+xi1)*(xi-beps))
        c0 = (A+xi1)*c0*X/((B+xi1)*xi)&
          - ((A-1.0)*(xn+2.*xi-1.0)+xi*(xi-beps))&
          *b0/(xi*(B+xi1)*(A+xi1-beps))
        t = c0 + xeps1*b0
        CHU = CHU + t
        IF ( ABS(t)<eps*ABS(CHU) ) RETURN
      END DO
      CALL XERMSG('SLATEC','CHU',&
        'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES',3,2)
    END IF
    !
    ! X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
    ! FORMULATION IS STABLE.
    !
    a0 = factor*pochai*GAMR(B+xi)*gamri1/beps
    b0 = xtoeps*b0/beps
    !
    CHU = summ + a0 - b0
    DO i = 1, 1000
      xi = istrt + i
      xi1 = istrt + i - 1
      a0 = (A+xi1)*a0*X/((B+xi1)*xi)
      b0 = (A+xi1-beps)*b0*X/((aintb+xi1)*(xi-beps))
      t = a0 - b0
      CHU = CHU + t
      IF ( ABS(t)<eps*ABS(CHU) ) RETURN
    END DO
    CALL XERMSG('SLATEC','CHU',&
      'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES',3,2)
  END IF
  !
  ! USE LUKE-S RATIONAL APPROX IN THE ASYMPTOTIC REGION.
  !
  CHU = X**(-A)*R9CHU(A,B,X)
  !
  RETURN
END FUNCTION CHU
