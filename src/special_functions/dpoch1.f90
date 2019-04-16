!** DPOCH1
REAL(8) FUNCTION DPOCH1(A,X)
  !>
  !***
  !  Calculate a generalization of Pochhammer's symbol starting
  !            from first order.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C1, C7A
  !***
  ! **Type:**      DOUBLE PRECISION (POCH1-S, DPOCH1-D)
  !***
  ! **Keywords:**  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate a double precision generalization of Pochhammer's symbol
  ! for double precision A and X for special situations that require
  ! especially accurate values when X is small in
  !        POCH1(A,X) = (POCH(A,X)-1)/X
  !                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
  ! This specification is particularly suited for stably computing
  ! expressions such as
  !        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
  !             = POCH1(A,X) - POCH1(B,X)
  ! Note that POCH1(A,0.0) = PSI(A)
  !
  ! When ABS(X) is so small that substantial cancellation will occur if
  ! the straightforward formula is used, we use an expansion due
  ! to Fields and discussed by Y. L. Luke, The Special Functions and Their
  ! Approximations, Vol. 1, Academic Press, 1969, page 34.
  !
  ! The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
  !        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
  ! In order to maintain significance in POCH1, we write for positive a
  !        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
  !                       = 1.0 + Q*EXPREL(Q) .
  ! Likewise the polynomial is written
  !        POLY = 1.0 + X*POLY1(A,X) .
  ! Thus,
  !        POCH1(A,X) = (POCH(A,X) - 1) / X
  !                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCOT, DEXPRL, DPOCH, DPSI, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)

  INTEGER i, ii, incr, j, k, ndx, nterms
  REAL(8) :: A, X, absa, absx, alnvar, b, binv, bp, gbern(21), &
    gbk, poly1, q, rho, sinpxx, sinpx2, term, trig, var, var2
  REAL(8), SAVE :: sqtbig, alneps
  REAL(8), PARAMETER :: bern(20) = [ +.833333333333333333333333333333333D-1, &
    -.138888888888888888888888888888888D-2, +.330687830687830687830687830687830D-4, &
    -.826719576719576719576719576719576D-6, +.208767569878680989792100903212014D-7, &
    -.528419013868749318484768220217955D-9, +.133825365306846788328269809751291D-10, &
    -.338968029632258286683019539124944D-12, +.858606205627784456413590545042562D-14, &
    -.217486869855806187304151642386591D-15, +.550900282836022951520265260890225D-17, &
    -.139544646858125233407076862640635D-18, +.353470703962946747169322997780379D-20, &
    -.895351742703754685040261131811274D-22, +.226795245233768306031095073886816D-23, &
    -.574472439520264523834847971943400D-24, +.145517247561486490186626486727132D-26, &
    -.368599494066531017818178247990866D-28, +.933673425709504467203255515278562D-30, &
    -.236502241570062993455963519636983D-31 ]
  REAL(8), PARAMETER :: pi = 3.141592653589793238462643383279503D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DPOCH1
  IF ( first ) THEN
    sqtbig = 1.0D0/SQRT(24.0D0*D1MACH(1))
    alneps = LOG(D1MACH(3))
    first = .FALSE.
  END IF
  !
  IF ( X==0.0D0 ) THEN
    DPOCH1 = DPSI(A)
    RETURN
  END IF
  !
  absx = ABS(X)
  absa = ABS(A)
  IF ( absx>0.1D0*absa ) THEN
    !
    DPOCH1 = (DPOCH(A,X)-1.0D0)/X
  ELSEIF ( absx*LOG(MAX(absa,2.0D0))>0.1D0 ) THEN
    DPOCH1 = (DPOCH(A,X)-1.0D0)/X
  ELSE
    !
    bp = A
    IF ( A<(-0.5D0) ) bp = 1.0D0 - A - X
    incr = 0
    IF ( bp<10.0D0 ) incr = 11 - INT( bp )
    b = bp + incr
    !
    var = b + 0.5D0*(X-1.0D0)
    alnvar = LOG(var)
    q = X*alnvar
    !
    poly1 = 0.0D0
    IF ( var<sqtbig ) THEN
      var2 = (1.0D0/var)**2
      !
      rho = 0.5D0*(X+1.0D0)
      gbern(1) = 1.0D0
      gbern(2) = -rho/12.0D0
      term = var2
      poly1 = gbern(2)*term
      !
      nterms = INT( -0.5D0*alneps/alnvar ) + 1
      IF ( nterms>20 ) CALL XERMSG('SLATEC','DPOCH1',&
        'NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD',1,2)
      IF ( nterms>=2 ) THEN
        !
        DO k = 2, nterms
          gbk = 0.0D0
          DO j = 1, k
            ndx = k - j + 1
            gbk = gbk + bern(ndx)*gbern(j)
          END DO
          gbern(k+1) = -rho*gbk/k
          !
          term = term*(2*k-2-X)*(2*k-1-X)*var2
          poly1 = poly1 + gbern(k+1)*term
        END DO
      END IF
    END IF
    !
    poly1 = (X-1.0D0)*poly1
    DPOCH1 = DEXPRL(q)*(alnvar+q*poly1) + poly1
    !
    IF ( incr/=0 ) THEN
      !
      ! WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
      ! TO OBTAIN DPOCH1(BP,X).
      !
      DO ii = 1, incr
        i = incr - ii
        binv = 1.0D0/(bp+i)
        DPOCH1 = (DPOCH1-binv)/(1.0D0+X*binv)
      END DO
    END IF
    !
    IF ( bp==A ) RETURN
    !
    ! WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
    ! FORMULA TO OBTAIN DPOCH1(A,X).
    !
    sinpxx = SIN(pi*X)/X
    sinpx2 = SIN(0.5D0*pi*X)
    trig = sinpxx*DCOT(pi*b) - 2.0D0*sinpx2*(sinpx2/X)
    !
    DPOCH1 = trig + (1.0D0+X*trig)*DPOCH1
    RETURN
  END IF
  !
END FUNCTION DPOCH1
