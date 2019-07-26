!** DEXINT
PURE SUBROUTINE DEXINT(X,N,Kode,M,Tol,En,Nz,Ierr)
  !> Compute an M member sequence of exponential integrals
  !     E(N+K,X), K=0,1,...,M-1 for N >= 1 and X >= 0.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      DOUBLE PRECISION (EXINT-S, DEXINT-D)
  !***
  ! **Keywords:**  EXPONENTIAL INTEGRAL, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !         DEXINT computes M member sequences of exponential integrals
  !         E(N+K,X), K=0,1,...,M-1 for N >= 1 and X >= 0.  The
  !         exponential integral is defined by
  !
  !         E(N,X)=integral on (1,infinity) of EXP(-XT)/T**N
  !
  !         where X=0.0 and N=1 cannot occur simultaneously.  Formulas
  !         and notation are found in the NBS Handbook of Mathematical
  !         Functions (ref. 1).
  !
  !         The power series is implemented for X <= XCUT and the
  !         confluent hypergeometric representation
  !
  !                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)
  !
  !         is computed for X > XCUT.  Since sequences are computed in
  !         a stable fashion by recurring away from X, A is selected as
  !         the integer closest to X within the constraint N <= A <=
  !         N+M-1.  For the U computation, A is further modified to be the
  !         nearest even integer.  Indices are carried forward or
  !         backward by the two term recursion relation
  !
  !                     K*E(K+1,X) + X*E(K,X) = EXP(-X)
  !
  !         once E(A,X) is computed.  The U function is computed by means
  !         of the backward recursive Miller algorithm applied to the
  !         three term contiguous relation for U(A+K,A,X), K=0,1,...
  !         This produces accurate ratios and determines U(A+K,A,X), and
  !         hence E(A,X), to within a multiplicative constant C.
  !         Another contiguous relation applied to C*U(A,A,X) and
  !         C*U(A+1,A,X) gets C*U(A+1,A+1,X), a quantity proportional to
  !         E(A+1,X).  The normalizing constant C is obtained from the
  !         two term recursion relation above with K=A.
  !
  !         The maximum number of significant digits obtainable
  !         is the smaller of 14 and the number of digits carried in
  !         double precision arithmetic.
  !
  !     Description of Arguments
  !
  !         Input     * X and TOL are double precision *
  !           X       X > 0.0 for N=1 and  X >= 0.0 for N >= 2
  !           N       order of the first member of the sequence, N >= 1
  !                   (X=0.0 and N=1 is an error)
  !           KODE    a selection parameter for scaled values
  !                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1.
  !                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1.
  !           M       number of exponential integrals in the sequence,
  !                   M >= 1
  !           TOL     relative accuracy wanted, ETOL <= TOL <= 0.1
  !                   ETOL is the larger of double precision unit
  !                   roundoff = eps_dp and 1.0D-18
  !
  !         Output    * EN is a double precision vector *
  !           EN      a vector of dimension at least M containing values
  !                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M
  !                   depending on KODE
  !           NZ      underflow indicator
  !                   NZ=0   a normal return
  !                   NZ=M   X exceeds XLIM and an underflow occurs.
  !                          EN(K)=0.0D0, K=1,M returned on KODE=1
  !           IERR    error flag
  !                   IERR=0, normal return, computation completed
  !                   IERR=1, input error,   no computation
  !                   IERR=2, error,         no computation
  !                           algorithm termination condition not met
  !
  !***
  ! **References:**  M. Abramowitz and I. A. Stegun, Handbook of
  !                 Mathematical Functions, NBS AMS Series 55, U.S. Dept.
  !                 of Commerce, 1955.
  !               D. E. Amos, Computation of exponential integrals, ACM
  !                 Transactions on Mathematical Software 6, (1980),
  !                 pp. 365-377 and pp. 420-428.
  !***
  ! **Routines called:**  D1MACH, DPSIXN, I1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   910408  Updated the REFERENCES section.  (WRB)
  !   920207  Updated with code with a revision date of 880811 from
  !           D. Amos.  Included correction of argument list.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : eps_dp, log10_radix_dp, min_exp_dp
  !
  INTEGER, INTENT(IN) :: Kode, M, N
  INTEGER, INTENT(OUT) :: Ierr, Nz
  REAL(DP), INTENT(IN) :: X, Tol
  REAL(DP), INTENT(OUT) :: En(M)
  !
  REAL(DP) :: a(99), aa, aams, ah, ak, at, b(99), bk, bt, cc, cnorm, ct, em, emx, &
    etol, fnm, fx, pt, p1, p2, s, tx, xlim, xtol, y(2), yt, y1, y2
  INTEGER :: i, ic, icase, ict, ik, ind, ix, i1m, jset, k, kk, kn, ks, ml, mu, nd, nm
  REAL(DP), PARAMETER :: xcut = 2._DP
  !* FIRST EXECUTABLE STATEMENT  DEXINT
  Ierr = 0
  Nz = 0
  etol = MAX(eps_dp,0.5E-18_DP)
  IF( X<0._DP ) Ierr = 1
  IF( N<1 ) Ierr = 1
  IF( Kode<1 .OR. Kode>2 ) Ierr = 1
  IF( M<1 ) Ierr = 1
  IF( Tol<etol .OR. Tol>0.1_DP ) Ierr = 1
  IF( X==0._DP .AND. N==1 ) Ierr = 1
  IF( Ierr/=0 ) RETURN
  i1m = -min_exp_dp
  pt = 2.3026_DP*i1m*log10_radix_dp
  xlim = pt - 6.907755_DP
  bt = pt + (N+M-1)
  IF( bt>1000._DP ) xlim = pt - LOG(bt)
  !
  IF( X>xcut ) THEN
    !-----------------------------------------------------------------------
    !     BACKWARD RECURSIVE MILLER ALGORITHM FOR
    !              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
    !     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
    !     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
    !-----------------------------------------------------------------------
    emx = 1._DP
    IF( Kode/=2 ) THEN
      IF( X<=xlim ) THEN
        emx = EXP(-X)
      ELSE
        Nz = M
        DO i = 1, M
          En(i) = 0._DP
        END DO
        RETURN
      END IF
    END IF
    tx = X + 0.5_DP
    ix = INT( tx )
    kn = N + M - 1
    IF( kn<=ix ) THEN
      icase = 1
      ks = kn
      ml = M - 1
      mu = -1
      ind = M
      IF( kn>1 ) GOTO 200
    ELSE
      IF( N<ix .AND. ix<kn ) GOTO 100
      IF( N>=ix ) THEN
        icase = 2
        ind = 1
        ks = N
        mu = M - 1
        IF( N>1 ) GOTO 200
        IF( kn/=1 ) THEN
          ix = 2
          GOTO 100
        END IF
      ELSE
        Ierr = 2
        RETURN
      END IF
    END IF
    ks = 2
    icase = 3
    GOTO 200
  ELSEIF( X==0._DP .AND. N>1 ) THEN
    DO i = 1, M
      En(i) = 1._DP/(N+i-2)
    END DO
    RETURN
  ELSE
    !-----------------------------------------------------------------------
    !     SERIES FOR E(N,X) FOR X<=XCUT
    !-----------------------------------------------------------------------
    tx = X + 0.5_DP
    ix = INT( tx )
    !-----------------------------------------------------------------------
    !     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
    !     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N>=2
    !-----------------------------------------------------------------------
    icase = 2
    IF( ix>N ) icase = 1
    nm = N - icase + 1
    nd = nm + 1
    ind = 3 - icase
    mu = M - ind
    ml = 1
    ks = nd
    fnm = nm
    s = 0._DP
    xtol = 3._DP*Tol
    IF( nd/=1 ) THEN
      xtol = 0.3333_DP*Tol
      s = 1._DP/fnm
    END IF
    aa = 1._DP
    ak = 1._DP
    ic = 35
    IF( X<etol ) ic = 1
    DO i = 1, ic
      aa = -aa*X/ak
      IF( i==nm ) THEN
        s = s + aa*(-LOG(X)+DPSIXN(nd))
        xtol = 3._DP*Tol
      ELSE
        s = s - aa/(ak-fnm)
        IF( ABS(aa)>xtol*ABS(s) ) THEN
          ak = ak + 1._DP
          CYCLE
        ELSEIF( i>=2 ) THEN
          IF( nd-2>i .OR. i>nd-1 ) GOTO 50
          ak = ak + 1._DP
          CYCLE
        END IF
      END IF
      ak = ak + 1._DP
    END DO
    IF( ic/=1 ) THEN
      Ierr = 2
      RETURN
    END IF
    50  IF( nd==1 ) s = s + (-LOG(X)+DPSIXN(1))
    IF( Kode==2 ) s = s*EXP(X)
    En(1) = s
    emx = 1._DP
    IF( M/=1 ) THEN
      En(ind) = s
      aa = ks
      IF( Kode==1 ) emx = EXP(-X)
      IF( icase==1 ) GOTO 300
      IF( icase==2 ) GOTO 400
    END IF
    IF( icase==2 ) RETURN
    IF( Kode==1 ) emx = EXP(-X)
    En(1) = (emx-s)/X
    RETURN
  END IF
  100  icase = 1
  ks = ix
  ml = ix - N
  ind = ml + 1
  mu = kn - ix
  200  ik = ks/2
  ah = ik
  jset = 1 + ks - (ik+ik)
  !-----------------------------------------------------------------------
  !     START COMPUTATION FOR
  !              EN(IND) = C*U( A, A ,X)    JSET=1
  !              EN(IND) = C*U(A+1,A+1,X)    JSET=2
  !     FOR AN EVEN INTEGER A.
  !-----------------------------------------------------------------------
  ic = 0
  aa = ah + ah
  aams = aa - 1._DP
  aams = aams*aams
  tx = X + X
  fx = tx + tx
  ak = ah
  xtol = Tol
  IF( Tol<=1.0D-3 ) xtol = 20._DP*Tol
  ct = aams + fx*ah
  em = (ah+1._DP)/((X+aa)*xtol*SQRT(ct))
  bk = aa
  cc = ah*ah
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
  !     RECURSION
  !-----------------------------------------------------------------------
  p1 = 0._DP
  p2 = 1._DP
  DO WHILE( ic/=99 )
    ic = ic + 1
    ak = ak + 1._DP
    at = bk/(bk+ak+cc+ic)
    bk = bk + ak + ak
    a(ic) = at
    bt = (ak+ak+X)/(ak+1._DP)
    b(ic) = bt
    pt = p2
    p2 = bt*p2 - at*p1
    p1 = pt
    ct = ct + fx
    em = em*at*(1._DP-tx/ct)
    IF( em*(ak+1._DP)<=p1*p1 ) THEN
      ict = ic
      kk = ic + 1
      bt = tx/(ct+fx)
      y2 = (bk/(bk+cc+kk))*(p1/p2)*(1._DP-bt+0.375_DP*bt*bt)
      y1 = 1._DP
      !-----------------------------------------------------------------------
      !     BACKWARD RECURRENCE FOR
      !              Y1=             C*U( A ,A,X)
      !              Y2= C*(A/(1+A/2))*U(A+1,A,X)
      !-----------------------------------------------------------------------
      DO k = 1, ict
        kk = kk - 1
        yt = y1
        y1 = (b(kk)*y1-y2)/a(kk)
        y2 = yt
      END DO
      !-----------------------------------------------------------------------
      !     THE CONTIGUOUS RELATION
      !              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
      !     WITH  B=A+1, C=A IS USED FOR
      !              Y(2) = C * U(A+1,A+1,X)
      !     X IS INCORPORATED INTO THE NORMALIZING RELATION
      !-----------------------------------------------------------------------
      pt = y2/y1
      cnorm = 1._SP - pt*(ah+1._SP)/aa
      y(1) = 1._SP/(cnorm*aa+X)
      y(2) = cnorm*y(1)
      IF( icase/=3 ) THEN
        En(ind) = emx*y(jset)
        IF( M==1 ) RETURN
        aa = ks
        IF( icase==1 ) GOTO 300
        IF( icase==2 ) GOTO 400
      END IF
      !-----------------------------------------------------------------------
      !     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX
      !-----------------------------------------------------------------------
      En(1) = emx*(1._SP-y(1))/X
      RETURN
    END IF
  END DO
  Ierr = 2
  RETURN
  300  k = ind - 1
  DO i = 1, ml
    aa = aa - 1._DP
    En(k) = (emx-aa*En(k+1))/X
    k = k - 1
  END DO
  IF( mu<=0 ) RETURN
  aa = ks
  400  k = ind
  DO i = 1, mu
    En(k+1) = (emx-X*En(k))/aa
    aa = aa + 1._DP
    k = k + 1
  END DO

  RETURN
END SUBROUTINE DEXINT