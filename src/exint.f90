!*==EXINT.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK EXINT
SUBROUTINE EXINT(X,N,Kode,M,Tol,En,Nz,Ierr)
  IMPLICIT NONE
  !*--EXINT5
  !***BEGIN PROLOGUE  EXINT
  !***PURPOSE  Compute an M member sequence of exponential integrals
  !            E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.
  !***LIBRARY   SLATEC
  !***CATEGORY  C5
  !***TYPE      SINGLE PRECISION (EXINT-S, DEXINT-D)
  !***KEYWORDS  EXPONENTIAL INTEGRAL, SPECIAL FUNCTIONS
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !         EXINT computes M member sequences of exponential integrals
  !         E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.  The
  !         exponential integral is defined by
  !
  !         E(N,X)=integral on (1,infinity) of EXP(-XT)/T**N
  !
  !         where X=0.0 and N=1 cannot occur simultaneously.  Formulas
  !         and notation are found in the NBS Handbook of Mathematical
  !         Functions (ref. 1).
  !
  !         The power series is implemented for X .LE. XCUT and the
  !         confluent hypergeometric representation
  !
  !                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)
  !
  !         is computed for X .GT. XCUT.  Since sequences are computed in
  !         a stable fashion by recurring away from X, A is selected as
  !         the integer closest to X within the constraint N .LE. A .LE.
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
  !     Description of Arguments
  !
  !         Input
  !           X       X .GT. 0.0 for N=1 and  X .GE. 0.0 for N .GE. 2
  !           N       order of the first member of the sequence, N .GE. 1
  !                   (X=0.0 and N=1 is an error)
  !           KODE    a selection parameter for scaled values
  !                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1.
  !                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1.
  !           M       number of exponential integrals in the sequence,
  !                   M .GE. 1
  !           TOL     relative accuracy wanted, ETOL .LE. TOL .LE. 0.1
  !                   ETOL = single precision unit roundoff = R1MACH(4)
  !
  !         Output
  !           EN      a vector of dimension at least M containing values
  !                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M
  !                   depending on KODE
  !           NZ      underflow indicator
  !                   NZ=0   a normal return
  !                   NZ=M   X exceeds XLIM and an underflow occurs.
  !                          EN(K)=0.0E0, K=1,M returned on KODE=1
  !           IERR    error flag
  !                   IERR=0, normal return, computation completed
  !                   IERR=1, input error,   no computation
  !                   IERR=2, error,         no computation
  !                           algorithm termination condition not met
  !
  !***REFERENCES  M. Abramowitz and I. A. Stegun, Handbook of
  !                 Mathematical Functions, NBS AMS Series 55, U.S. Dept.
  !                 of Commerce, 1955.
  !               D. E. Amos, Computation of exponential integrals, ACM
  !                 Transactions on Mathematical Software 6, (1980),
  !                 pp. 365-377 and pp. 420-428.
  !***ROUTINES CALLED  I1MACH, PSIXN, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   910408  Updated the REFERENCES section.  (WRB)
  !   920207  Updated with code with a revision date of 880811 from
  !           D. Amos.  Included correction of argument list.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  EXINT
  REAL a, aa, aams, ah, ak, at, b, bk, bt, cc, cnorm, ct, em, &
    emx, En, etol, fnm, fx, pt, p1, p2, s, Tol, tx, X, xcut, &
    xlim, xtol, y, yt, y1, y2
  REAL R1MACH, PSIXN
  INTEGER i, ic, icase, ict, Ierr, ik, ind, ix, i1m, jset, k, &
    kk, kn, Kode, ks, M, ml, mu, N, nd, nm, Nz
  INTEGER I1MACH
  DIMENSION En(*), a(99), b(99), y(2)
  !***FIRST EXECUTABLE STATEMENT  EXINT
  Ierr = 0
  Nz = 0
  etol = MAX(R1MACH(4),0.5E-18)
  IF ( X<0.0E0 ) Ierr = 1
  IF ( N<1 ) Ierr = 1
  IF ( Kode<1.OR.Kode>2 ) Ierr = 1
  IF ( M<1 ) Ierr = 1
  IF ( Tol<etol.OR.Tol>0.1E0 ) Ierr = 1
  IF ( X==0.0E0.AND.N==1 ) Ierr = 1
  IF ( Ierr/=0 ) RETURN
  i1m = -I1MACH(12)
  pt = 2.3026E0*R1MACH(5)*i1m
  xlim = pt - 6.907755E0
  bt = pt + (N+M-1)
  IF ( bt>1000.0E0 ) xlim = pt - LOG(bt)
  !
  xcut = 2.0E0
  IF ( etol>2.0E-7 ) xcut = 1.0E0
  IF ( X>xcut ) THEN
    !-----------------------------------------------------------------------
    !     BACKWARD RECURSIVE MILLER ALGORITHM FOR
    !              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
    !     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
    !     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
    !-----------------------------------------------------------------------
    emx = 1.0E0
    IF ( Kode/=2 ) THEN
      IF ( X<=xlim ) THEN
        emx = EXP(-X)
      ELSE
        Nz = M
        DO i = 1, M
          En(i) = 0.0E0
        ENDDO
        RETURN
      ENDIF
    ENDIF
    ix = X + 0.5E0
    kn = N + M - 1
    IF ( kn<=ix ) THEN
      icase = 1
      ks = kn
      ml = M - 1
      mu = -1
      ind = M
      IF ( kn>1 ) GOTO 200
    ELSE
      IF ( N<ix.AND.ix<kn ) GOTO 100
      IF ( N>=ix ) THEN
        icase = 2
        ind = 1
        ks = N
        mu = M - 1
        IF ( N>1 ) GOTO 200
        IF ( kn/=1 ) THEN
          ix = 2
          GOTO 100
        ENDIF
      ELSE
        Ierr = 2
        GOTO 99999
      ENDIF
    ENDIF
    ks = 2
    icase = 3
    GOTO 200
  ELSEIF ( X==0.0E0.AND.N>1 ) THEN
    DO i = 1, M
      En(i) = 1.0E0/(N+i-2)
    ENDDO
    RETURN
  ELSE
    !-----------------------------------------------------------------------
    !     SERIES FOR E(N,X) FOR X.LE.XCUT
    !-----------------------------------------------------------------------
    tx = X + 0.5E0
    ix = tx
    !-----------------------------------------------------------------------
    !     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
    !     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N.GE.2
    !-----------------------------------------------------------------------
    icase = 2
    IF ( ix>N ) icase = 1
    nm = N - icase + 1
    nd = nm + 1
    ind = 3 - icase
    mu = M - ind
    ml = 1
    ks = nd
    fnm = nm
    s = 0.0E0
    xtol = 3.0E0*Tol
    IF ( nd/=1 ) THEN
      xtol = 0.3333E0*Tol
      s = 1.0E0/fnm
    ENDIF
    aa = 1.0E0
    ak = 1.0E0
    ic = 35
    IF ( X<etol ) ic = 1
    DO i = 1, ic
      aa = -aa*X/ak
      IF ( i==nm ) THEN
        s = s + aa*(-LOG(X)+PSIXN(nd))
        xtol = 3.0E0*Tol
      ELSE
        s = s - aa/(ak-fnm)
        IF ( ABS(aa)>xtol*ABS(s) ) THEN
          ak = ak + 1.0E0
          CYCLE
        ELSEIF ( i>=2 ) THEN
          IF ( nd-2>i.OR.i>nd-1 ) GOTO 50
          ak = ak + 1.0E0
          CYCLE
        ENDIF
      ENDIF
      ak = ak + 1.0E0
    ENDDO
    IF ( ic/=1 ) THEN
      Ierr = 2
      GOTO 99999
    ENDIF
    50     IF ( nd==1 ) s = s + (-LOG(X)+PSIXN(1))
    IF ( Kode==2 ) s = s*EXP(X)
    En(1) = s
    emx = 1.0E0
    IF ( M/=1 ) THEN
      En(ind) = s
      aa = ks
      IF ( Kode==1 ) emx = EXP(-X)
      IF ( icase==1 ) GOTO 300
      IF ( icase==2 ) GOTO 400
    ENDIF
    IF ( icase==2 ) RETURN
    IF ( Kode==1 ) emx = EXP(-X)
    En(1) = (emx-s)/X
    RETURN
  ENDIF
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
  aams = aa - 1.0E0
  aams = aams*aams
  tx = X + X
  fx = tx + tx
  ak = ah
  xtol = Tol
  IF ( Tol<=1.0E-3 ) xtol = 20.0E0*Tol
  ct = aams + fx*ah
  em = (ah+1.0E0)/((X+aa)*xtol*SQRT(ct))
  bk = aa
  cc = ah*ah
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
  !     RECURSION
  !-----------------------------------------------------------------------
  p1 = 0.0E0
  p2 = 1.0E0
  DO WHILE ( ic/=99 )
    ic = ic + 1
    ak = ak + 1.0E0
    at = bk/(bk+ak+cc+ic)
    bk = bk + ak + ak
    a(ic) = at
    bt = (ak+ak+X)/(ak+1.0E0)
    b(ic) = bt
    pt = p2
    p2 = bt*p2 - at*p1
    p1 = pt
    ct = ct + fx
    em = em*at*(1.0E0-tx/ct)
    IF ( em*(ak+1.0E0)<=p1*p1 ) THEN
      ict = ic
      kk = ic + 1
      bt = tx/(ct+fx)
      y2 = (bk/(bk+cc+kk))*(p1/p2)*(1.0E0-bt+0.375E0*bt*bt)
      y1 = 1.0E0
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
      ENDDO
      !-----------------------------------------------------------------------
      !     THE CONTIGUOUS RELATION
      !              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
      !     WITH  B=A+1, C=A IS USED FOR
      !              Y(2) = C * U(A+1,A+1,X)
      !     X IS INCORPORATED INTO THE NORMALIZING RELATION
      !-----------------------------------------------------------------------
      pt = y2/y1
      cnorm = 1.0E0 - pt*(ah+1.0E0)/aa
      y(1) = 1.0E0/(cnorm*aa+X)
      y(2) = cnorm*y(1)
      IF ( icase/=3 ) THEN
        En(ind) = emx*y(jset)
        IF ( M==1 ) RETURN
        aa = ks
        IF ( icase==1 ) GOTO 300
        IF ( icase==2 ) GOTO 400
      ENDIF
      !-----------------------------------------------------------------------
      !     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX
      !-----------------------------------------------------------------------
      En(1) = emx*(1.0E0-y(1))/X
      RETURN
    ENDIF
  ENDDO
  Ierr = 2
  GOTO 99999
  300  k = ind - 1
  DO i = 1, ml
    aa = aa - 1.0E0
    En(k) = (emx-aa*En(k+1))/X
    k = k - 1
  ENDDO
  IF ( mu<=0 ) RETURN
  aa = ks
  400  k = ind
  DO i = 1, mu
    En(k+1) = (emx-X*En(k))/aa
    aa = aa + 1.0E0
    k = k + 1
  ENDDO
  RETURN
  99999 CONTINUE
  END SUBROUTINE EXINT
