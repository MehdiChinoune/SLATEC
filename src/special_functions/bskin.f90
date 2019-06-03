!** BSKIN
SUBROUTINE BSKIN(X,N,Kode,M,Y,Nz,Ierr)
  !>
  !  Compute repeated integrals of the K-zero Bessel function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10F
  !***
  ! **Type:**      SINGLE PRECISION (BSKIN-S, DBSKIN-D)
  !***
  ! **Keywords:**  BICKLEY FUNCTIONS, EXPONENTIAL INTEGRAL,
  !             INTEGRALS OF BESSEL FUNCTIONS, K-ZERO BESSEL FUNCTION
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !         The following definitions are used in BSKIN:
  !
  !   Definition 1
  !         KI(0,X) = K-zero Bessel function.
  !
  !   Definition 2
  !         KI(N,X) = Bickley Function
  !                 =  integral from X to infinity of KI(N-1,t)dt
  !                     for X .ge. 0 and N = 1,2,...
  !   ____________________________________________________________________
  !      BSKIN computes sequences of Bickley functions (repeated integrals
  !      of the K0 Bessel function); i.e. for fixed X and N and K=1,...,
  !      BSKIN computes the M-member sequence
  !
  !                     Y(K) =        KI(N+K-1,X) for KODE=1
  !      or
  !                     Y(K) = EXP(X)*KI(N+K-1,X) for KODE=2,
  !
  !      for N.ge.0 and X.ge.0 (N and X cannot be zero simultaneously).
  !
  !      INPUT
  !        X      - Argument, X .ge. 0.0E0
  !        N      - Order of first member of the sequence N .ge. 0
  !        KODE   - Selection parameter
  !                 KODE = 1 returns Y(K)=       KI(N+K-1,X), K=1,M
  !                      = 2 returns Y(K)=EXP(X)*KI(N+K-1,X), K=1,M
  !        M      - Number of members in the sequence, M.ge.1
  !
  !      OUTPUT
  !        Y      - A vector of dimension at least M containing the
  !                 sequence selected by KODE.
  !        NZ     - Underflow flag
  !                 NZ = 0 means computation completed
  !                    = M means an exponential underflow occurred on
  !                        KODE=1.  Y(K)=0.0E0, K=1,...,M is returned
  !        IERR   - Error flag
  !                 IERR = 0, Normal return, computation completed.
  !                      = 1, Input error,   no computation.
  !                      = 2, Error,         no computation.  The
  !                           termination condition was not met.
  !
  !      The nominal computational accuracy is the maximum of unit
  !      roundoff (=R1MACH(4)) and 1.0e-18 since critical constants
  !      are given to only 18 digits.
  !
  !      DBSKIN is the double precision version of BSKIN.
  !
  !- Long Description:
  !
  !         Numerical recurrence on
  !
  !      (L-1)*KI(L,X) = X(KI(L-3,X) - KI(L-1,X)) + (L-2)*KI(L-2,X)
  !
  !         is stable where recurrence is carried forward or backward
  !         away from INT(X+0.5).  The power series for indices 0,1 and 2
  !         on 0.le.X.le. 2 starts a stable recurrence for indices
  !         greater than 2.  If N is sufficiently large (N.gt.NLIM), the
  !         uniform asymptotic expansion for N to INFINITY is more
  !         economical.  On X.gt.2 the recursion is started by evaluating
  !         the uniform expansion for the three members whose indices are
  !         closest to INT(X+0.5) within the set N,...,N+M-1.  Forward
  !         recurrence, backward recurrence or both, complete the
  !         sequence depending on the relation of INT(X+0.5) to the
  !         indices N,...,N+M-1.
  !
  !***
  ! **References:**  D. E. Amos, Uniform asymptotic expansions for
  !                 exponential integrals E(N,X) and Bickley functions
  !                 KI(N,X), ACM Transactions on Mathematical Software,
  !                 1983.
  !               D. E. Amos, A portable Fortran subroutine for the
  !                 Bickley functions KI(N,X), Algorithm 609, ACM
  !                 Transactions on Mathematical Software, 1983.
  !***
  ! **Routines called:**  BKIAS, BKISR, EXINT, GAMRN, I1MACH, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891009  Removed unreferenced statement label.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH, I1MACH
  INTEGER :: Ierr, Kode, M, N, Nz
  REAL(SP) :: X, Y(M)
  INTEGER :: i, icase, il, i1m, k, kk, ktrms, m3, ne, nflg, nl, nlim, nn, np, ns, nt
  REAL(SP) enlim, exi(102), fn, gr, h(31), hn, ss, tol, t1, t2, w, &
    xlim, xnlim, xp, ys(3), yss(3)
  !-----------------------------------------------------------------------
  !             COEFFICIENTS IN SERIES OF EXPONENTIAL INTEGRALS
  !-----------------------------------------------------------------------
  REAL(SP), PARAMETER :: a(50) = [ 1.00000000000000000E+00, 5.00000000000000000E-01, &
    3.75000000000000000E-01, 3.12500000000000000E-01, 2.73437500000000000E-01, &
    2.46093750000000000E-01, 2.25585937500000000E-01, 2.09472656250000000E-01, &
    1.96380615234375000E-01, 1.85470581054687500E-01, 1.76197052001953125E-01, &
    1.68188095092773438E-01, 1.61180257797241211E-01, 1.54981017112731934E-01, &
    1.49445980787277222E-01, 1.44464448094367981E-01, 1.39949934091418982E-01, &
    1.35833759559318423E-01, 1.32060599571559578E-01, 1.28585320635465905E-01, &
    1.25370687619579257E-01, 1.22385671247684513E-01, 1.19604178719328047E-01, &
    1.17004087877603524E-01, 1.14566502713486784E-01, 1.12275172659217048E-01, &
    1.10116034723462874E-01, 1.08076848895250599E-01, 1.06146905164978267E-01, &
    1.04316786110409676E-01, 1.02578173008569515E-01, 1.00923686347140974E-01, &
    9.93467537479668965E-02, 9.78414999033007314E-02, 9.64026543164874854E-02, &
    9.50254735405376642E-02, 9.37056752969190855E-02, 9.24393823875012600E-02, &
    9.12230747245078224E-02, 9.00535481254756708E-02, 8.89278787739072249E-02, &
    8.78433924473961612E-02, 8.67976377754033498E-02, 8.57883629175498224E-02, &
    8.48134951571231199E-02, 8.38711229887106408E-02, 8.29594803475290034E-02, &
    8.20769326842574183E-02, 8.12219646354630702E-02, 8.03931690779583449E-02 ]
  !-----------------------------------------------------------------------
  !             SQRT(PI)/2
  !-----------------------------------------------------------------------
  REAL(SP), PARAMETER :: hrtpi = 8.86226925452758014E-01
  !
  !* FIRST EXECUTABLE STATEMENT  BSKIN
  Ierr = 0
  Nz = 0
  IF ( X<0.0E0 ) Ierr = 1
  IF ( N<0 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( M<1 ) Ierr = 1
  IF ( X==0.0E0.AND.N==0 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  IF ( X==0.0E0 ) THEN
    !-----------------------------------------------------------------------
    !     X=0 CASE
    !-----------------------------------------------------------------------
    fn = N
    hn = 0.5E0*fn
    gr = GAMRN(hn)
    Y(1) = hrtpi*gr
    IF ( M==1 ) RETURN
    Y(2) = hrtpi/(hn*gr)
    IF ( M==2 ) RETURN
    DO k = 3, M
      Y(k) = fn*Y(k-2)/(fn+1.0E0)
      fn = fn + 1.0E0
    END DO
    RETURN
  ELSE
    i1m = -I1MACH(12)
    t1 = 2.3026E0*R1MACH(5)*i1m
    xlim = t1 - 3.228086E0
    t2 = t1 + N + M - 1
    IF ( t2>1000.0E0 ) xlim = t1 - 0.5E0*(LOG(t2)-0.451583E0)
    IF ( X>xlim.AND.Kode==1 ) GOTO 400
    tol = MAX(R1MACH(4),1.0E-18)
    i1m = I1MACH(11)
    !-----------------------------------------------------------------------
    !     LN(NLIM) = 0.125*LN(EPS),   NLIM = 2*KTRMS+N
    !-----------------------------------------------------------------------
    xnlim = 0.287823E0*(i1m-1)*R1MACH(5)
    enlim = EXP(xnlim)
    nlim = INT(enlim) + 2
    nlim = MIN(100,nlim)
    nlim = MAX(20,nlim)
    m3 = MIN(M,3)
    nl = N + M - 1
    IF ( X>2.0E0 ) THEN
      !-----------------------------------------------------------------------
      !     COMPUTATION BY ASYMPTOTIC EXPANSION FOR X.GT.2
      !-----------------------------------------------------------------------
      w = X + 0.5E0
      nt = INT(w)
      IF ( nl>nt ) THEN
        IF ( N>=nt ) GOTO 300
        !-----------------------------------------------------------------------
        !     ICASE=2, N.LT.NT.LT.NL WITH BOTH FORWARD AND BACKWARD RECURSION
        !-----------------------------------------------------------------------
        nn = nt + 1
        nflg = MIN(M-m3,1)
        icase = 2
        GOTO 200
      ELSE
        !-----------------------------------------------------------------------
        !     CASE NL.LE.NT, ICASE=0
        !-----------------------------------------------------------------------
        icase = 0
        nn = nl
        nflg = MIN(M-m3,1)
        GOTO 200
      END IF
    ELSE
      IF ( N>nlim ) GOTO 300
      !-----------------------------------------------------------------------
      !     COMPUTATION BY SERIES FOR 0.LE.X.LE.2
      !-----------------------------------------------------------------------
      nflg = 0
      nn = N
      IF ( nl>2 ) THEN
        m3 = 3
        nn = 0
        nflg = 1
      END IF
      xp = 1.0E0
      IF ( Kode==2 ) xp = EXP(X)
      DO i = 1, m3
        CALL BKISR(X,nn,w,Ierr)
        IF ( Ierr/=0 ) RETURN
        w = w*xp
        IF ( nn>=N ) THEN
          kk = nn - N + 1
          Y(kk) = w
        END IF
        ys(i) = w
        nn = nn + 1
      END DO
      IF ( nflg==0 ) RETURN
      ns = nn
      xp = 1.0E0
    END IF
  END IF
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION SCALED BY EXP(X) ON ICASE=0,1,2
  !-----------------------------------------------------------------------
  100  fn = ns - 1
  il = nl - ns + 1
  IF ( il<=0 ) RETURN
  DO i = 1, il
    t1 = ys(2)
    t2 = ys(3)
    ys(3) = (X*(ys(1)-ys(3))+(fn-1.0E0)*ys(2))/fn
    ys(2) = t2
    ys(1) = t1
    fn = fn + 1.0E0
    IF ( ns>=N ) THEN
      kk = ns - N + 1
      Y(kk) = ys(3)*xp
    END IF
    ns = ns + 1
  END DO
  RETURN
  200  kk = (nlim-nn)/2
  ktrms = MAX(0,kk)
  ns = nn + 1
  np = nn - m3 + 1
  xp = 1.0E0
  IF ( Kode==1 ) xp = EXP(-X)
  DO i = 1, m3
    kk = i
    CALL BKIAS(X,np,ktrms,a,w,kk,ne,gr,h,Ierr)
    IF ( Ierr/=0 ) RETURN
    ys(i) = w
    np = np + 1
  END DO
  !-----------------------------------------------------------------------
  !     SUM SERIES OF EXPONENTIAL INTEGRALS BACKWARD
  !-----------------------------------------------------------------------
  IF ( ktrms/=0 ) THEN
    ne = ktrms + ktrms + 1
    np = nn - m3 + 2
    CALL EXINT(X,np,2,ne,tol,exi,Nz,Ierr)
    IF ( Nz/=0 ) GOTO 400
    IF ( Ierr==2 ) RETURN
  END IF
  DO i = 1, m3
    ss = 0.0E0
    IF ( ktrms/=0 ) THEN
      kk = i + ktrms + ktrms - 2
      il = ktrms
      DO k = 1, ktrms
        ss = ss + a(il)*exi(kk)
        kk = kk - 2
        il = il - 1
      END DO
    END IF
    ys(i) = ys(i) + ss
  END DO
  IF ( icase/=1 ) THEN
    IF ( nflg/=0 ) THEN
      !-----------------------------------------------------------------------
      !     BACKWARD RECURSION SCALED BY EXP(X) ICASE=0,2
      !-----------------------------------------------------------------------
      kk = nn - N + 1
      k = m3
      DO i = 1, m3
        Y(kk) = ys(k)*xp
        yss(i) = ys(i)
        kk = kk - 1
        k = k - 1
      END DO
      il = kk
      IF ( il>0 ) THEN
        fn = nn - 3
        DO i = 1, il
          t1 = ys(2)
          t2 = ys(1)
          ys(1) = ys(2) + ((fn+2.0E0)*ys(3)-(fn+1.0E0)*ys(1))/X
          ys(2) = t2
          ys(3) = t1
          Y(kk) = ys(1)*xp
          kk = kk - 1
          fn = fn - 1.0E0
        END DO
      END IF
      IF ( icase/=2 ) RETURN
      DO i = 1, m3
        ys(i) = yss(i)
      END DO
      GOTO 100
    END IF
  END IF
  DO i = 1, m3
    Y(i) = ys(i)*xp
  END DO
  IF ( icase==1.AND.nflg==1 ) GOTO 100
  RETURN
  !-----------------------------------------------------------------------
  !     ICASE=1, NT.LE.N.LE.NL WITH FORWARD RECURSION
  !-----------------------------------------------------------------------
  300  nn = N + m3 - 1
  nflg = MIN(M-m3,1)
  icase = 1
  GOTO 200
  !-----------------------------------------------------------------------
  !     UNDERFLOW ON KODE=1, X.GT.XLIM
  !-----------------------------------------------------------------------
  400  Nz = M
  DO i = 1, M
    Y(i) = 0.0E0
  END DO
END SUBROUTINE BSKIN
