!** DPOLFT
SUBROUTINE DPOLFT(N,X,Y,W,Maxdeg,Ndeg,Eps,R,Ierr,A)
  !>
  !  Fit discrete data in a least squares sense by polynomials
  !            in one variable.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1A1A2
  !***
  ! **Type:**      DOUBLE PRECISION (POLFIT-S, DPOLFT-D)
  !***
  ! **Keywords:**  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
  !***
  ! **Author:**  Shampine, L. F., (SNLA)
  !           Davenport, S. M., (SNLA)
  !           Huddleston, R. E., (SNLL)
  !***
  ! **Description:**
  !
  !     Abstract
  !
  !     Given a collection of points X(I) and a set of values Y(I) which
  !     correspond to some function or measurement at each of the X(I),
  !     subroutine  DPOLFT  computes the weighted least-squares polynomial
  !     fits of all degrees up to some degree either specified by the user
  !     or determined by the routine.  The fits thus obtained are in
  !     orthogonal polynomial form.  Subroutine  DP1VLU  may then be
  !     called to evaluate the fitted polynomials and any of their
  !     derivatives at any point.  The subroutine  DPCOEF  may be used to
  !     express the polynomial fits as powers of (X-C) for any specified
  !     point C.
  !
  !     The parameters for  DPOLFT  are
  !
  !     Input -- All TYPE REAL variables are DOUBLE PRECISION
  !         N -      the number of data points.  The arrays X, Y and W
  !                  must be dimensioned at least  N  (N .GE. 1).
  !         X -      array of values of the independent variable.  These
  !                  values may appear in any order and need not all be
  !                  distinct.
  !         Y -      array of corresponding function values.
  !         W -      array of positive values to be used as weights.  If
  !                  W(1) is negative,  DPOLFT  will set all the weights
  !                  to 1.0, which means unweighted least squares error
  !                  will be minimized.  To minimize relative error, the
  !                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
  !                  I = 1,...,N .
  !         MAXDEG - maximum degree to be allowed for polynomial fit.
  !                  MAXDEG  may be any non-negative integer less than  N.
  !                  Note -- MAXDEG  cannot be equal to  N-1  when a
  !                  statistical test is to be used for degree selection,
  !                  i.e., when input value of  EPS  is negative.
  !         EPS -    specifies the criterion to be used in determining
  !                  the degree of fit to be computed.
  !                  (1)  If  EPS  is input negative,  DPOLFT  chooses the
  !                       degree based on a statistical F test of
  !                       significance.  One of three possible
  !                       significance levels will be used:  .01, .05 or
  !                       .10.  If  EPS=-1.0, the routine will
  !                       automatically select one of these levels based
  !                       on the number of data points and the maximum
  !                       degree to be considered.  If  EPS  is input as
  !                       -.01, -.05, or -.10, a significance level of
  !                       .01, .05, or .10, respectively, will be used.
  !                  (2)  If  EPS  is set to 0.,  DPOLFT  computes the
  !                       polynomials of degrees 0 through  MAXDEG .
  !                  (3)  If  EPS  is input positive,  EPS  is the RMS
  !                       error tolerance which must be satisfied by the
  !                       fitted polynomial.  DPOLFT  will increase the
  !                       degree of fit until this criterion is met or
  !                       until the maximum degree is reached.
  !
  !     Output -- All TYPE REAL variables are DOUBLE PRECISION
  !         NDEG -   degree of the highest degree fit computed.
  !         EPS -    RMS error of the polynomial of degree  NDEG .
  !         R -      vector of dimension at least NDEG containing values
  !                  of the fit of degree  NDEG  at each of the  X(I) .
  !                  Except when the statistical test is used, these
  !                  values are more accurate than results from subroutine
  !                  DP1VLU  normally are.
  !         IERR -   error flag with the following possible values.
  !             1 -- indicates normal execution, i.e., either
  !                  (1)  the input value of  EPS  was negative, and the
  !                       computed polynomial fit of degree  NDEG
  !                       satisfies the specified F test, or
  !                  (2)  the input value of  EPS  was 0., and the fits of
  !                       all degrees up to  MAXDEG  are complete, or
  !                  (3)  the input value of  EPS  was positive, and the
  !                       polynomial of degree  NDEG  satisfies the RMS
  !                       error requirement.
  !             2 -- invalid input parameter.  At least one of the input
  !                  parameters has an illegal value and must be corrected
  !                  before  DPOLFT  can proceed.  Valid input results
  !                  when the following restrictions are observed
  !                       N .GE. 1
  !                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
  !                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
  !                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
  !             3 -- cannot satisfy the RMS error requirement with a
  !                  polynomial of degree no greater than  MAXDEG .  Best
  !                  fit found is of degree  MAXDEG .
  !             4 -- cannot satisfy the test for significance using
  !                  current value of  MAXDEG .  Statistically, the
  !                  best fit found is of order  NORD .  (In this case,
  !                  NDEG will have one of the values:  MAXDEG-2,
  !                  MAXDEG-1, or MAXDEG).  Using a higher value of
  !                  MAXDEG  may result in passing the test.
  !         A -      work and output array having at least 3N+3MAXDEG+3
  !                  locations
  !
  !     Note - DPOLFT  calculates all fits of degrees up to and including
  !            NDEG .  Any or all of these fits can be evaluated or
  !            expressed as powers of (X-C) using  DP1VLU  and  DPCOEF
  !            after just one call to  DPOLFT .
  !
  !***
  ! **References:**  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
  !                 Curve fitting by polynomials in one variable, Report
  !                 SLA-74-0270, Sandia Laboratories, June 1974.
  !***
  ! **Routines called:**  DP1VLU, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   740601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900911  Added variable YP to DOUBLE PRECISION declaration.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
  USE service, ONLY : XERMSG
  INTEGER N, i, idegf, Ierr, j, jp1, jpas, k1, k1pj, k2, k2pj, k3, &
    k3pi, k4, k4pi, k5, k5pi, ksig, m, Maxdeg, mop1, Ndeg, nder, nfail
  REAL(8) :: temd1, temd2, A(*), degf, den, Eps, etst, f, fcrit, R(*), sig, &
    sigj, sigjm1, sigpas, temp, X(*), xm, Y(*), yp(1), W(*), w1, w11
  REAL(8), PARAMETER :: co(4,3)= RESHAPE( [ &
    -13.086850D0, -2.4648165D0, -3.3846535D0, -1.2973162D0, &
    -3.3381146D0, -1.7812271D0, -3.2578406D0, -1.6589279D0, &
    -1.6282703D0, -1.3152745D0, -3.2640179D0, -1.9829776D0 ], [4,3] )
  !* FIRST EXECUTABLE STATEMENT  DPOLFT
  m = ABS(N)
  yp = 0.D0
  IF ( m==0 ) GOTO 700
  IF ( Maxdeg<0 ) GOTO 700
  A(1) = Maxdeg
  mop1 = Maxdeg + 1
  IF ( m<mop1 ) GOTO 700
  IF ( Eps<0.0D0.AND.m==mop1 ) GOTO 700
  xm = m
  etst = Eps*Eps*xm
  IF ( W(1)<0.0D0 ) THEN
    DO i = 1, m
      W(i) = 1.0D0
    END DO
  ELSE
    DO i = 1, m
      IF ( W(i)<=0.0D0 ) GOTO 700
    END DO
  END IF
  IF ( Eps<0.0D0 ) THEN
    !
    ! DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
    ! CHOOSING DEGREE OF POLYNOMIAL FIT
    !
    IF ( Eps>(-.55D0) ) THEN
      ksig = 1
      IF ( Eps<(-.03D0) ) ksig = 2
      IF ( Eps<(-.07D0) ) ksig = 3
    ELSE
      idegf = m - Maxdeg - 1
      ksig = 1
      IF ( idegf<10 ) ksig = 2
      IF ( idegf<5 ) ksig = 3
    END IF
  END IF
  !
  ! INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
  !
  k1 = Maxdeg + 1
  k2 = k1 + Maxdeg
  k3 = k2 + Maxdeg + 2
  k4 = k3 + m
  k5 = k4 + m
  DO i = 2, k4
    A(i) = 0.0D0
  END DO
  w11 = 0.0D0
  IF ( N<0 ) THEN
    !
    ! CONSTRAINED CASE
    !
    DO i = 1, m
      k4pi = k4 + i
      w11 = w11 + W(i)*A(k4pi)**2
    END DO
  ELSE
    !
    ! UNCONSTRAINED CASE
    !
    DO i = 1, m
      k4pi = k4 + i
      A(k4pi) = 1.0D0
      w11 = w11 + W(i)
    END DO
  END IF
  !
  ! COMPUTE FIT OF DEGREE ZERO
  !
  temd1 = 0.0D0
  DO i = 1, m
    k4pi = k4 + i
    temd1 = temd1 + W(i)*Y(i)*A(k4pi)
  END DO
  temd1 = temd1/w11
  A(k2+1) = temd1
  sigj = 0.0D0
  DO i = 1, m
    k4pi = k4 + i
    k5pi = k5 + i
    temd2 = temd1*A(k4pi)
    R(i) = temd2
    A(k5pi) = temd2 - R(i)
    sigj = sigj + W(i)*((Y(i)-R(i))-A(k5pi))**2
  END DO
  j = 0
  !
  ! SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
  !
  IF ( Eps<0 ) GOTO 200
  IF ( Eps==0 ) GOTO 300
  GOTO 400
  100 CONTINUE
  DO
    !
    ! INCREMENT DEGREE
    !
    j = j + 1
    jp1 = j + 1
    k1pj = k1 + j
    k2pj = k2 + j
    sigjm1 = sigj
    !
    ! COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
    !
    IF ( j>1 ) A(k1pj) = w11/w1
    !
    ! COMPUTE NEW A COEFFICIENT
    !
    temd1 = 0.0D0
    DO i = 1, m
      k4pi = k4 + i
      temd2 = A(k4pi)
      temd1 = temd1 + X(i)*W(i)*temd2*temd2
    END DO
    A(jp1) = temd1/w11
    !
    ! EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
    !
    w1 = w11
    w11 = 0.0D0
    DO i = 1, m
      k3pi = k3 + i
      k4pi = k4 + i
      temp = A(k3pi)
      A(k3pi) = A(k4pi)
      A(k4pi) = (X(i)-A(jp1))*A(k3pi) - A(k1pj)*temp
      w11 = w11 + W(i)*A(k4pi)**2
    END DO
    !
    ! GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
    ! PRECISION
    !
    temd1 = 0.0D0
    DO i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = W(i)*((Y(i)-R(i))-A(k5pi))*A(k4pi)
      temd1 = temd1 + temd2
    END DO
    temd1 = temd1/w11
    A(k2pj+1) = temd1
    !
    ! UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
    ! ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
    ! COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
    ! THE MOST SIGNIFICANT BITS ARE STORED IN  R(I), AND THE LEAST
    ! SIGNIFICANT BITS ARE IN  A(K5PI) .
    !
    sigj = 0.0D0
    DO i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = R(i) + A(k5pi) + temd1*A(k4pi)
      R(i) = temd2
      A(k5pi) = temd2 - R(i)
      sigj = sigj + W(i)*((Y(i)-R(i))-A(k5pi))**2
    END DO
    !
    ! SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
    ! MAXDEG  HAS BEEN REACHED
    !
    IF ( Eps<0 ) THEN
      !
      ! COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
      !
      IF ( sigj==0.0D0 ) GOTO 600
      degf = m - j - 1
      den = (co(4,ksig)*degf+1.0D0)*degf
      fcrit = (((co(3,ksig)*degf)+co(2,ksig))*degf+co(1,ksig))/den
      fcrit = fcrit*fcrit
      f = (sigjm1-sigj)*degf/sigj
      IF ( f>=fcrit ) EXIT
      !
      ! POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
      ! SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
      !
      nfail = nfail + 1
      IF ( nfail>=3 ) GOTO 600
      IF ( Maxdeg==j ) GOTO 800
    ELSEIF ( Eps==0 ) THEN
      GOTO 300
    ELSE
      GOTO 400
    END IF
  END DO
  !
  ! POLYNOMIAL OF DEGREE J SATISFIES F TEST
  !
  200  sigpas = sigj
  jpas = j
  nfail = 0
  IF ( Maxdeg/=j ) GOTO 100
  GOTO 800
  !
  ! RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
  ! EPS = 0.)
  !
  300 CONTINUE
  IF ( Maxdeg/=j ) GOTO 100
  GOTO 500
  !
  ! SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
  !
  400 CONTINUE
  IF ( sigj>etst ) THEN
    IF ( Maxdeg/=j ) GOTO 100
    Ierr = 3
    Ndeg = Maxdeg
    sig = sigj
    GOTO 900
  END IF
  !
  ! RETURNS
  !
  500  Ierr = 1
  Ndeg = j
  sig = sigj
  GOTO 900
  600  Ierr = 1
  Ndeg = jpas
  sig = sigpas
  GOTO 900
  700  Ierr = 2
  CALL XERMSG('SLATEC','DPOLFT','INVALID INPUT PARAMETER.',2,1)
  RETURN
  800  Ierr = 4
  Ndeg = jpas
  sig = sigpas
  !
  900  A(k3) = Ndeg
  !
  ! WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
  ! ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
  !
  IF ( Eps<0.0.AND.Ndeg/=Maxdeg ) THEN
    nder = 0
    DO i = 1, m
      CALL DP1VLU(Ndeg,nder,X(i),R(i),yp,A)
    END DO
  END IF
  Eps = SQRT(sig/xm)
  RETURN
END SUBROUTINE DPOLFT
