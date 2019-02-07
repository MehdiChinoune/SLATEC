!*==PSIFN.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PSIFN
SUBROUTINE PSIFN(X,N,Kode,M,Ans,Nz,Ierr)
  IMPLICIT NONE
  !*--PSIFN5
  !***BEGIN PROLOGUE  PSIFN
  !***PURPOSE  Compute derivatives of the Psi function.
  !***LIBRARY   SLATEC
  !***CATEGORY  C7C
  !***TYPE      SINGLE PRECISION (PSIFN-S, DPSIFN-D)
  !***KEYWORDS  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
  !             PSI FUNCTION
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !         The following definitions are used in PSIFN:
  !
  !      Definition 1
  !         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
  !                  the LOG GAMMA function.
  !      Definition 2
  !                     K   K
  !         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
  !   ___________________________________________________________________
  !       PSIFN computes a sequence of SCALED derivatives of
  !       the PSI function; i.e. for fixed X and M it computes
  !       the M-member sequence
  !
  !                  ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
  !                    for K = N,...,N+M-1
  !
  !       where PSI(K,X) is as defined above.   For KODE=1, PSIFN returns
  !       the scaled derivatives as described.  KODE=2 is operative only
  !       when K=0 and in that case PSIFN returns -PSI(X) + LN(X).  That
  !       is, the logarithmic behavior for large X is removed when KODE=1
  !       and K=0.  When sums or differences of PSI functions are computed
  !       the logarithmic terms can be combined analytically and computed
  !       separately to help retain significant digits.
  !
  !         Note that CALL PSIFN(X,0,1,1,ANS) results in
  !                   ANS = -PSI(X)
  !
  !     Input
  !           X      - Argument, X .gt. 0.0E0
  !           N      - First member of the sequence, 0 .le. N .le. 100
  !                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
  !                                       -PSI(X)+LN(X) for KODE=2
  !           KODE   - Selection parameter
  !                    KODE=1 returns scaled derivatives of the PSI
  !                    function.
  !                    KODE=2 returns scaled derivatives of the PSI
  !                    function EXCEPT when N=0. In this case,
  !                    ANS(1) = -PSI(X) + LN(X) is returned.
  !           M      - Number of members of the sequence, M .ge. 1
  !
  !    Output
  !           ANS    - A vector of length at least M whose first M
  !                    components contain the sequence of derivatives
  !                    scaled according to KODE.
  !           NZ     - Underflow flag
  !                    NZ.eq.0, A normal return
  !                    NZ.ne.0, Underflow, last NZ components of ANS are
  !                             set to zero, ANS(M-K+1)=0.0, K=1,...,NZ
  !           IERR   - Error flag
  !                    IERR=0, A normal return, computation completed
  !                    IERR=1, Input error,     no computation
  !                    IERR=2, Overflow,        X too small or N+M-1 too
  !                            large or both
  !                    IERR=3, Error,           N too large. Dimensioned
  !                            array TRMR(NMAX) is not large enough for N
  !
  !         The nominal computational accuracy is the maximum of unit
  !         roundoff (=R1MACH(4)) and 1.0E-18 since critical constants
  !         are given to only 18 digits.
  !
  !         DPSIFN is the Double Precision version of PSIFN.
  !
  ! *Long Description:
  !
  !         The basic method of evaluation is the asymptotic expansion
  !         for large X.ge.XMIN followed by backward recursion on a two
  !         term recursion relation
  !
  !                  W(X+1) + X**(-N-1) = W(X).
  !
  !         This is supplemented by a series
  !
  !                  SUM( (X+K)**(-N-1) , K=0,1,2,... )
  !
  !         which converges rapidly for large N. Both XMIN and the
  !         number of terms of the series are calculated from the unit
  !         roundoff of the machine environment.
  !
  !***REFERENCES  Handbook of Mathematical Functions, National Bureau
  !                 of Standards Applied Mathematics Series 55, edited
  !                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
  !                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
  !               D. E. Amos, A portable Fortran subroutine for
  !                 derivatives of the Psi function, Algorithm 610, ACM
  !                 Transactions on Mathematical Software 9, 4 (1983),
  !                 pp. 494-502.
  !***ROUTINES CALLED  I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  PSIFN
  INTEGER i , Ierr , j , k , Kode , M , mm , mx , N , nmax , nn , np , nx , &
    Nz
  INTEGER I1MACH
  REAL Ans , arg , b , den , elim , eps , fln , fn , fnp , fns , fx , rln , &
    rxsq , r1m4 , r1m5 , s , slope , t , ta , tk , tol , tols , trm , &
    trmr , tss , tst , tt , t1 , t2 , wdtol , X , xdmln , xdmy , xinc , &
    xln , xm , xmin , xq , yint
  REAL R1MACH
  DIMENSION b(22) , trm(22) , trmr(100) , Ans(*)
  SAVE nmax , b
  DATA nmax/100/
  !-----------------------------------------------------------------------
  !             BERNOULLI NUMBERS
  !-----------------------------------------------------------------------
  DATA b(1) , b(2) , b(3) , b(4) , b(5) , b(6) , b(7) , b(8) , b(9) , &
    b(10) , b(11) , b(12) , b(13) , b(14) , b(15) , b(16) , b(17) , &
    b(18) , b(19) , b(20) , b(21) , b(22)/1.00000000000000000E+00 , &
    -5.00000000000000000E-01 , 1.66666666666666667E-01 , &
    -3.33333333333333333E-02 , 2.38095238095238095E-02 , &
    -3.33333333333333333E-02 , 7.57575757575757576E-02 , &
    -2.53113553113553114E-01 , 1.16666666666666667E+00 , &
    -7.09215686274509804E+00 , 5.49711779448621554E+01 , &
    -5.29124242424242424E+02 , 6.19212318840579710E+03 , &
    -8.65802531135531136E+04 , 1.42551716666666667E+06 , &
    -2.72982310678160920E+07 , 6.01580873900642368E+08 , &
    -1.51163157670921569E+10 , 4.29614643061166667E+11 , &
    -1.37116552050883328E+13 , 4.88332318973593167E+14 , &
    -1.92965793419400681E+16/
  !
  !***FIRST EXECUTABLE STATEMENT  PSIFN
  Ierr = 0
  Nz = 0
  IF ( X<=0.0E0 ) Ierr = 1
  IF ( N<0 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( M<1 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  mm = M
  nx = MIN(-I1MACH(12),I1MACH(13))
  r1m5 = R1MACH(5)
  r1m4 = R1MACH(4)*0.5E0
  wdtol = MAX(r1m4,0.5E-18)
  !-----------------------------------------------------------------------
  !     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  elim = 2.302E0*(nx*r1m5-3.0E0)
  xln = LOG(X)
  100  nn = N + mm - 1
  fn = nn
  fnp = fn + 1.0E0
  t = fnp*xln
  !-----------------------------------------------------------------------
  !     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X
  !-----------------------------------------------------------------------
  IF ( ABS(t)<=elim ) THEN
    IF ( X<wdtol ) THEN
      !-----------------------------------------------------------------------
      !     SMALL X.LT.UNIT ROUND OFF
      !-----------------------------------------------------------------------
      Ans(1) = X**(-N-1)
      IF ( mm/=1 ) THEN
        k = 1
        DO i = 2 , mm
          Ans(k+1) = Ans(k)/X
          k = k + 1
        ENDDO
      ENDIF
      IF ( N/=0 ) RETURN
      IF ( Kode==2 ) Ans(1) = Ans(1) + xln
      RETURN
    ELSE
      !-----------------------------------------------------------------------
      !     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1
      !-----------------------------------------------------------------------
      rln = r1m5*I1MACH(11)
      rln = MIN(rln,18.06E0)
      fln = MAX(rln,3.0E0) - 3.0E0
      yint = 3.50E0 + 0.40E0*fln
      slope = 0.21E0 + fln*(0.0006038E0*fln+0.008677E0)
      xm = yint + slope*fn
      mx = INT(xm) + 1
      xmin = mx
      IF ( N/=0 ) THEN
        xm = -2.302E0*rln - MIN(0.0E0,xln)
        fns = N
        arg = xm/fns
        arg = MIN(0.0E0,arg)
        eps = EXP(arg)
        xm = 1.0E0 - eps
        IF ( ABS(arg)<1.0E-3 ) xm = -arg
        fln = X*xm/eps
        xm = xmin - X
        IF ( xm>7.0E0.AND.fln<15.0E0 ) THEN
          !-----------------------------------------------------------------------
          !     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,...
          !-----------------------------------------------------------------------
          nn = INT(fln) + 1
          np = N + 1
          t1 = (fns+1.0E0)*xln
          t = EXP(-t1)
          s = t
          den = X
          DO i = 1 , nn
            den = den + 1.0E0
            trm(i) = den**(-np)
            s = s + trm(i)
          ENDDO
          Ans(1) = s
          IF ( N==0 ) THEN
            IF ( Kode==2 ) Ans(1) = s + xln
          ENDIF
          IF ( mm==1 ) RETURN
          !-----------------------------------------------------------------------
          !     GENERATE HIGHER DERIVATIVES, J.GT.N
          !-----------------------------------------------------------------------
          tol = wdtol/5.0E0
          DO j = 2 , mm
            t = t/X
            s = t
            tols = t*tol
            den = X
            DO i = 1 , nn
              den = den + 1.0E0
              trm(i) = trm(i)/den
              s = s + trm(i)
              IF ( trm(i)<tols ) EXIT
            ENDDO
            Ans(j) = s
          ENDDO
          RETURN
        ENDIF
      ENDIF
      xdmy = X
      xdmln = xln
      xinc = 0.0E0
      IF ( X<xmin ) THEN
        nx = INT(X)
        xinc = xmin - nx
        xdmy = X + xinc
        xdmln = LOG(xdmy)
      ENDIF
      !-----------------------------------------------------------------------
      !     GENERATE W(N+MM-1,X) BY THE ASYMPTOTIC EXPANSION
      !-----------------------------------------------------------------------
      t = fn*xdmln
      t1 = xdmln + xdmln
      t2 = t + xdmln
      tk = MAX(ABS(t),ABS(t1),ABS(t2))
      IF ( tk>elim ) GOTO 200
      tss = EXP(-t)
      tt = 0.5E0/xdmy
      t1 = tt
      tst = wdtol*tt
      IF ( nn/=0 ) t1 = tt + 1.0E0/fn
      rxsq = 1.0E0/(xdmy*xdmy)
      ta = 0.5E0*rxsq
      t = fnp*ta
      s = t*b(3)
      IF ( ABS(s)>=tst ) THEN
        tk = 2.0E0
        DO k = 4 , 22
          t = t*((tk+fn+1.0E0)/(tk+1.0E0))*((tk+fn)/(tk+2.0E0))*rxsq
          trm(k) = t*b(k)
          IF ( ABS(trm(k))<tst ) EXIT
          s = s + trm(k)
          tk = tk + 2.0E0
        ENDDO
      ENDIF
      s = (s+t1)*tss
      IF ( xinc/=0.0E0 ) THEN
        !-----------------------------------------------------------------------
        !     BACKWARD RECUR FROM XDMY TO X
        !-----------------------------------------------------------------------
        nx = INT(xinc)
        np = nn + 1
        IF ( nx>nmax ) THEN
          Ierr = 3
          Nz = 0
          GOTO 99999
        ELSE
          IF ( nn==0 ) GOTO 120
          xm = xinc - 1.0E0
          fx = X + xm
          !-----------------------------------------------------------------------
          !     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL
          !-----------------------------------------------------------------------
          DO i = 1 , nx
            trmr(i) = fx**(-np)
            s = s + trmr(i)
            xm = xm - 1.0E0
            fx = X + xm
          ENDDO
        ENDIF
      ENDIF
      Ans(mm) = s
      IF ( fn==0.0E0 ) GOTO 150
      !-----------------------------------------------------------------------
      !     GENERATE LOWER DERIVATIVES, J.LT.N+MM-1
      !-----------------------------------------------------------------------
      IF ( mm==1 ) RETURN
      DO j = 2 , mm
        fnp = fn
        fn = fn - 1.0E0
        tss = tss*xdmy
        t1 = tt
        IF ( fn/=0.0E0 ) t1 = tt + 1.0E0/fn
        t = fnp*ta
        s = t*b(3)
        IF ( ABS(s)>=tst ) THEN
          tk = 3.0E0 + fnp
          DO k = 4 , 22
            trm(k) = trm(k)*fnp/tk
            IF ( ABS(trm(k))<tst ) EXIT
            s = s + trm(k)
            tk = tk + 2.0E0
          ENDDO
        ENDIF
        s = (s+t1)*tss
        IF ( xinc/=0.0E0 ) THEN
          IF ( fn==0.0E0 ) GOTO 120
          xm = xinc - 1.0E0
          fx = X + xm
          DO i = 1 , nx
            trmr(i) = trmr(i)*fx
            s = s + trmr(i)
            xm = xm - 1.0E0
            fx = X + xm
          ENDDO
        ENDIF
        mx = mm - j + 1
        Ans(mx) = s
        IF ( fn==0.0E0 ) GOTO 150
      ENDDO
      RETURN
      !-----------------------------------------------------------------------
      !     RECURSION FOR N = 0
      !-----------------------------------------------------------------------
      120      DO i = 1 , nx
      s = s + 1.0E0/(X+nx-i)
    ENDDO
    ENDIF
    150    IF ( Kode==2 ) THEN
    IF ( xdmy==X ) RETURN
    xq = xdmy/X
    Ans(1) = s - LOG(xq)
    RETURN
  ELSE
    Ans(1) = s - xdmln
    RETURN
  ENDIF
ELSEIF ( t<=0.0E0 ) THEN
  Nz = 0
  Ierr = 2
  RETURN
ENDIF
200  Nz = Nz + 1
Ans(mm) = 0.0E0
mm = mm - 1
IF ( mm==0 ) RETURN
GOTO 100
99999 END SUBROUTINE PSIFN
