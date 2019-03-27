!** DPSIFN
SUBROUTINE DPSIFN(X,N,Kode,M,Ans,Nz,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Compute derivatives of the Psi function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C7C
  !***
  ! **Type:**      DOUBLE PRECISION (PSIFN-S, DPSIFN-D)
  !***
  ! **Keywords:**  DERIVATIVES OF THE GAMMA FUNCTION, POLYGAMMA FUNCTION,
  !             PSI FUNCTION
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !         The following definitions are used in DPSIFN:
  !
  !      Definition 1
  !         PSI(X) = d/dx (ln(GAMMA(X)), the first derivative of
  !                  the log GAMMA function.
  !      Definition 2
  !                     K   K
  !         PSI(K,X) = d /dx (PSI(X)), the K-th derivative of PSI(X).
  !   ___________________________________________________________________
  !      DPSIFN computes a sequence of SCALED derivatives of
  !      the PSI function; i.e. for fixed X and M it computes
  !      the M-member sequence
  !
  !                    ((-1)**(K+1)/GAMMA(K+1))*PSI(K,X)
  !                       for K = N,...,N+M-1
  !
  !      where PSI(K,X) is as defined above.   For KODE=1, DPSIFN returns
  !      the scaled derivatives as described.  KODE=2 is operative only
  !      when K=0 and in that case DPSIFN returns -PSI(X) + LN(X).  That
  !      is, the logarithmic behavior for large X is removed when KODE=2
  !      and K=0.  When sums or differences of PSI functions are computed
  !      the logarithmic terms can be combined analytically and computed
  !      separately to help retain significant digits.
  !
  !         Note that CALL DPSIFN(X,0,1,1,ANS) results in
  !                   ANS = -PSI(X)
  !
  !     Input      X is DOUBLE PRECISION
  !           X      - Argument, X .gt. 0.0D0
  !           N      - First member of the sequence, 0 .le. N .le. 100
  !                    N=0 gives ANS(1) = -PSI(X)       for KODE=1
  !                                       -PSI(X)+LN(X) for KODE=2
  !           KODE   - Selection parameter
  !                    KODE=1 returns scaled derivatives of the PSI
  !                    function.
  !                    KODE=2 returns scaled derivatives of the PSI
  !                    function EXCEPT when N=0. In this case,
  !                    ANS(1) = -PSI(X) + LN(X) is returned.
  !           M      - Number of members of the sequence, M.ge.1
  !
  !    Output     ANS is DOUBLE PRECISION
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
  !         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants
  !         are given to only 18 digits.
  !
  !         PSIFN is the single precision version of DPSIFN.
  !
  !- Long Description:
  !
  !         The basic method of evaluation is the asymptotic expansion
  !         for large X.ge.XMIN followed by backward recursion on a two
  !         term recursion relation
  !
  !                  W(X+1) + X**(-N-1) = W(X).
  !
  !         This is supplemented by a series
  !
  !                  SUM( (X+K)**(-N-1), K=0,1,2,... )
  !
  !         which converges rapidly for large N. Both XMIN and the
  !         number of terms of the series are calculated from the unit
  !         roundoff of the machine environment.
  !
  !***
  ! **References:**  Handbook of Mathematical Functions, National Bureau
  !                 of Standards Applied Mathematics Series 55, edited
  !                 by M. Abramowitz and I. A. Stegun, equations 6.3.5,
  !                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
  !               D. E. Amos, A portable Fortran subroutine for
  !                 derivatives of the Psi function, Algorithm 610, ACM
  !                 Transactions on Mathematical Software 9, 4 (1983),
  !                 pp. 494-502.
  !***
  ! **Routines called:**  D1MACH, I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891006  Cosmetic changes to prologue.  (WRB)
  !   891006  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, Ierr, j, k, Kode, M, mm, mx, N, nmax, nn, np, nx, Nz, fn
  INTEGER I1MACH
  REAL(8) :: Ans(*), arg, den, elim, eps, fln, fx, rln, &
    rxsq, r1m4, r1m5, s, slope, t, ta, tk, tol, &
    tols, trm(22), trmr(100), tss, tst, tt, t1, t2, wdtol, &
    X, xdmln, xdmy, xinc, xln, xm, xmin, xq, yint
  REAL(8) :: D1MACH
  SAVE nmax
  DATA nmax/100/
  !-----------------------------------------------------------------------
  !             BERNOULLI NUMBERS
  !-----------------------------------------------------------------------
  REAL(8), PARAMETER :: b(22) = [ 1.00000000000000000D+00, -5.00000000000000000D-01, &
    1.66666666666666667D-01, -3.33333333333333333D-02, 2.38095238095238095D-02, &
    -3.33333333333333333D-02, 7.57575757575757576D-02, -2.53113553113553114D-01, &
    1.16666666666666667D+00, -7.09215686274509804D+00, 5.49711779448621554D+01, &
    -5.29124242424242424D+02, 6.19212318840579710D+03, -8.65802531135531136D+04, &
    1.42551716666666667D+06, -2.72982310678160920D+07, 6.01580873900642368D+08, &
    -1.51163157670921569D+10, 4.29614643061166667D+11, -1.37116552050883328D+13, &
    4.88332318973593167D+14, -1.92965793419400681D+16 ]
  !
  !* FIRST EXECUTABLE STATEMENT  DPSIFN
  Ierr = 0
  Nz = 0
  IF ( X<=0.0D0 ) Ierr = 1
  IF ( N<0 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( M<1 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  mm = M
  nx = MIN(-I1MACH(15),I1MACH(16))
  r1m5 = D1MACH(5)
  r1m4 = D1MACH(4)*0.5D0
  wdtol = MAX(r1m4,0.5D-18)
  !-----------------------------------------------------------------------
  !     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  elim = 2.302D0*(nx*r1m5-3.0D0)
  xln = LOG(X)
  100  nn = N + mm - 1
  fn = nn
  t = (fn+1)*xln
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
        DO i = 2, mm
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
      rln = r1m5*I1MACH(14)
      rln = MIN(rln,18.06D0)
      fln = MAX(rln,3.0D0) - 3.0D0
      yint = 3.50D0 + 0.40D0*fln
      slope = 0.21D0 + fln*(0.0006038D0*fln+0.008677D0)
      xm = yint + slope*fn
      mx = INT(xm) + 1
      xmin = mx
      IF ( N/=0 ) THEN
        xm = -2.302D0*rln - MIN(0.0D0,xln)
        arg = xm/N
        arg = MIN(0.0D0,arg)
        eps = EXP(arg)
        xm = 1.0D0 - eps
        IF ( ABS(arg)<1.0D-3 ) xm = -arg
        fln = X*xm/eps
        xm = xmin - X
        IF ( xm>7.0D0.AND.fln<15.0D0 ) THEN
          !-----------------------------------------------------------------------
          !     COMPUTE BY SERIES (X+K)**(-(N+1)), K=0,1,2,...
          !-----------------------------------------------------------------------
          nn = INT(fln) + 1
          np = N + 1
          t1 = (N+1)*xln
          t = EXP(-t1)
          s = t
          den = X
          DO i = 1, nn
            den = den + 1.0D0
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
          tol = wdtol/5.0D0
          DO j = 2, mm
            t = t/X
            s = t
            tols = t*tol
            den = X
            DO i = 1, nn
              den = den + 1.0D0
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
      xinc = 0.0D0
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
      tt = 0.5D0/xdmy
      t1 = tt
      tst = wdtol*tt
      IF ( nn/=0 ) t1 = tt + 1.0D0/fn
      rxsq = 1.0D0/(xdmy*xdmy)
      ta = 0.5D0*rxsq
      t = (fn+1)*ta
      s = t*b(3)
      IF ( ABS(s)>=tst ) THEN
        tk = 2.0D0
        DO k = 4, 22
          t = t*((tk+fn+1)/(tk+1.0D0))*((tk+fn)/(tk+2.0D0))*rxsq
          trm(k) = t*b(k)
          IF ( ABS(trm(k))<tst ) EXIT
          s = s + trm(k)
          tk = tk + 2.0D0
        ENDDO
      ENDIF
      s = (s+t1)*tss
      IF ( xinc/=0.0D0 ) THEN
        !-----------------------------------------------------------------------
        !     BACKWARD RECUR FROM XDMY TO X
        !-----------------------------------------------------------------------
        nx = INT(xinc)
        np = nn + 1
        IF ( nx>nmax ) THEN
          Nz = 0
          Ierr = 3
          RETURN
        ELSE
          IF ( nn==0 ) GOTO 120
          xm = xinc - 1.0D0
          fx = X + xm
          !-----------------------------------------------------------------------
          !     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL
          !-----------------------------------------------------------------------
          DO i = 1, nx
            trmr(i) = fx**(-np)
            s = s + trmr(i)
            xm = xm - 1.0D0
            fx = X + xm
          ENDDO
        ENDIF
      ENDIF
      Ans(mm) = s
      IF ( fn==0 ) GOTO 150
      !-----------------------------------------------------------------------
      !     GENERATE LOWER DERIVATIVES, J.LT.N+MM-1
      !-----------------------------------------------------------------------
      IF ( mm==1 ) RETURN
      DO j = 2, mm
        fn = fn - 1
        tss = tss*xdmy
        t1 = tt
        IF ( fn/=0 ) t1 = tt + 1.0D0/fn
        t = (fn+1)*ta
        s = t*b(3)
        IF ( ABS(s)>=tst ) THEN
          tk = 4 + fn
          DO k = 4, 22
            trm(k) = trm(k)*(fn+1)/tk
            IF ( ABS(trm(k))<tst ) EXIT
            s = s + trm(k)
            tk = tk + 2.0D0
          ENDDO
        ENDIF
        s = (s+t1)*tss
        IF ( xinc/=0.0D0 ) THEN
          IF ( fn==0 ) GOTO 120
          xm = xinc - 1.0D0
          fx = X + xm
          DO i = 1, nx
            trmr(i) = trmr(i)*fx
            s = s + trmr(i)
            xm = xm - 1.0D0
            fx = X + xm
          ENDDO
        ENDIF
        mx = mm - j + 1
        Ans(mx) = s
        IF ( fn==0 ) GOTO 150
      ENDDO
      RETURN
      !-----------------------------------------------------------------------
      !     RECURSION FOR N = 0
      !-----------------------------------------------------------------------
      120 CONTINUE
      DO i = 1, nx
        s = s + 1.0D0/(X+nx-i)
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
ELSEIF ( t<=0.0D0 ) THEN
  Nz = 0
  Ierr = 2
  RETURN
ENDIF
200  Nz = Nz + 1
Ans(mm) = 0.0D0
mm = mm - 1
IF ( mm==0 ) RETURN
GOTO 100
RETURN
END SUBROUTINE DPSIFN
