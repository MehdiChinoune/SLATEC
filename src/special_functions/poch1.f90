!** POCH1
REAL FUNCTION POCH1(A,X)
  IMPLICIT NONE
  !>
  !***
  !  Calculate a generalization of Pochhammer's symbol starting
  !            from first order.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C1, C7A
  !***
  ! **Type:**      SINGLE PRECISION (POCH1-S, DPOCH1-D)
  !***
  ! **Keywords:**  FIRST ORDER, FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate a generalization of Pochhammer's symbol for special
  ! situations that require especially accurate values when X is small in
  !        POCH1(A,X) = (POCH(A,X)-1)/X
  !                   = (GAMMA(A+X)/GAMMA(A) - 1.0)/X .
  ! This specification is particularly suited for stably computing
  ! expressions such as
  !        (GAMMA(A+X)/GAMMA(A) - GAMMA(B+X)/GAMMA(B))/X
  !             = POCH1(A,X) - POCH1(B,X)
  ! Note that POCH1(A,0.0) = PSI(A)
  !
  ! When ABS(X) is so small that substantial cancellation will occur if
  ! the straightforward formula is used, we  use an expansion due
  ! to Fields and discussed by Y. L. Luke, The Special Functions and Their
  ! Approximations, Vol. 1, Academic Press, 1969, page 34.
  !
  ! The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
  !        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
  ! In order to maintain significance in POCH1, we write for positive A
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
  ! **Routines called:**  COT, EXPREL, POCH, PSI, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)

  REAL A, absa, absx, alneps, alnvar, b, bern, binv, bp, gbern, gbk, pi, poly1, &
    q, rho, sinpx2, sinpxx, sqtbig, term, trig, var, var2, X
  INTEGER i, ii, incr, j, k, ndx, nterms
  DIMENSION bern(9), gbern(10)
  LOGICAL first
  REAL, EXTERNAL :: EXPREL, COT, POCH, PSI, R1MACH
  SAVE bern, pi, sqtbig, alneps, first
  DATA bern(1)/.83333333333333333E-01/
  DATA bern(2)/ - .13888888888888889E-02/
  DATA bern(3)/.33068783068783069E-04/
  DATA bern(4)/ - .82671957671957672E-06/
  DATA bern(5)/.20876756987868099E-07/
  DATA bern(6)/ - .52841901386874932E-09/
  DATA bern(7)/.13382536530684679E-10/
  DATA bern(8)/ - .33896802963225829E-12/
  DATA bern(9)/.85860620562778446E-14/
  DATA pi/3.14159265358979324E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  POCH1
  IF ( first ) THEN
    sqtbig = 1.0/SQRT(24.0*R1MACH(1))
    alneps = LOG(R1MACH(3))
  ENDIF
  first = .FALSE.
  !
  IF ( X==0.0 ) THEN
    POCH1 = PSI(A)
    RETURN
  ENDIF
  !
  absx = ABS(X)
  absa = ABS(A)
  IF ( absx>0.1*absa ) THEN
    !
    POCH1 = (POCH(A,X)-1.0)/X
  ELSEIF ( absx*LOG(MAX(absa,2.0))>0.1 ) THEN
    POCH1 = (POCH(A,X)-1.0)/X
  ELSE
    !
    bp = A
    IF ( A<(-0.5) ) bp = 1.0 - A - X
    incr = 0
    IF ( bp<10.0 ) incr = 11 - INT( bp )
    b = bp + incr
    !
    var = b + 0.5*(X-1.0)
    alnvar = LOG(var)
    q = X*alnvar
    !
    poly1 = 0.0
    IF ( var<sqtbig ) THEN
      var2 = (1.0/var)**2
      !
      rho = 0.5*(X+1.0)
      gbern(1) = 1.0
      gbern(2) = -rho/12.0
      term = var2
      poly1 = gbern(2)*term
      !
      nterms = INT( -0.5*alneps/alnvar ) + 1
      IF ( nterms>9 ) CALL XERMSG('SLATEC','POCH1',&
        'NTERMS IS TOO BIG, MAYBE R1MACH(3) IS BAD',1,2)
      IF ( nterms>=2 ) THEN
        !
        DO k = 2, nterms
          gbk = 0.0
          DO j = 1, k
            ndx = k - j + 1
            gbk = gbk + bern(ndx)*gbern(j)
          ENDDO
          gbern(k+1) = -rho*gbk/k
          !
          term = term*(2*k-2.-X)*(2*k-1.-X)*var2
          poly1 = poly1 + gbern(k+1)*term
        ENDDO
      ENDIF
    ENDIF
    !
    poly1 = (X-1.0)*poly1
    POCH1 = EXPREL(q)*(alnvar+q*poly1) + poly1
    !
    IF ( incr/=0 ) THEN
      !
      ! WE HAVE POCH1(B,X).  BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
      ! TO OBTAIN POCH1(BP,X).
      !
      DO ii = 1, incr
        i = incr - ii
        binv = 1.0/(bp+i)
        POCH1 = (POCH1-binv)/(1.0+X*binv)
      ENDDO
    ENDIF
    !
    IF ( bp==A ) RETURN
    !
    ! WE HAVE POCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A REFLECTION
    ! FORMULA TO OBTAIN POCH1(A,X).
    !
    sinpxx = SIN(pi*X)/X
    sinpx2 = SIN(0.5*pi*X)
    trig = sinpxx*COT(pi*b) - 2.0*sinpx2*(sinpx2/X)
    !
    POCH1 = trig + (1.0+X*trig)*POCH1
    RETURN
  ENDIF
  !
END FUNCTION POCH1
