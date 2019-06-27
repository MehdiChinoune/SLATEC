!** BESI
PURE SUBROUTINE BESI(X,Alpha,Kode,N,Y,Nz)
  !> Compute an N member sequence of I Bessel functions
  !  I_{ALPHA+K-1}(X), K=1,...,N or scaled Bessel functions
  !  EXP(-X)*I_{ALPHA+K-1}(X), K=1,...,N for non-negative ALPHA and X.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10B3
  !***
  ! **Type:**      SINGLE PRECISION (BESI-S, DBESI-D)
  !***
  ! **Keywords:**  I BESSEL FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !           Daniel, S. L., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BESI computes an N member sequence of I Bessel functions
  !         I_{ALPHA+K-1}(X), K=1,...,N or scaled Bessel functions
  !         EXP(-X)*I_{ALPHA+K-1}(X), K=1,...,N for non-negative ALPHA
  !         and X.  A combination of the power series, the asymptotic
  !         expansion for X to infinity, and the uniform asymptotic
  !         expansion for NU to infinity are applied over subdivisions of
  !         the (NU,X) plane.  For values not covered by one of these
  !         formulae, the order is incremented by an integer so that one
  !         of these formulae apply.  Backward recursion is used to reduce
  !         orders by integer values.  The asymptotic expansion for X to
  !         infinity is used only when the entire sequence (specifically
  !         the last member) lies within the region covered by the
  !         expansion.  Leading terms of these expansions are used to test
  !         for over or underflow where appropriate.  If a sequence is
  !         requested and the last member would underflow, the result is
  !         set to zero and the next lower order tried, etc., until a
  !         member comes on scale or all are set to zero.  An overflow
  !         cannot occur with scaling.
  !
  !     Description of Arguments
  !
  !         Input
  !           X      - X >= 0.0E0
  !           ALPHA  - order of first member of the sequence,
  !                    ALPHA >= 0.0E0
  !           KODE   - a parameter to indicate the scaling option
  !                    KODE=1 returns
  !                           Y(K)=        I_{ALPHA+K-1}(X),
  !                                K=1,...,N
  !                    KODE=2 returns
  !                           Y(K)=EXP(-X)*I_{ALPHA+K-1}(X),
  !                                K=1,...,N
  !           N      - number of members in the sequence, N >= 1
  !
  !         Output
  !           Y      - a vector whose first N components contain
  !                    values for I_{ALPHA+K-1}(X) or scaled
  !                    values for EXP(-X)*I_{ALPHA+K-1}(X),
  !                    K=1,...,N depending on KODE
  !           NZ     - number of components of Y set to zero due to
  !                    underflow,
  !                    NZ=0  , normal return, computation completed
  !                    NZ /= 0, last NZ components of Y set to zero,
  !                             Y(K)=0.0E0, K=N-NZ+1,...,N.
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Overflow with KODE=1 - a fatal error
  !         Underflow - a non-fatal error (NZ /= 0)
  !
  !***
  ! **References:**  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
  !                 subroutines IBESS and JBESS for Bessel functions
  !                 I(NU,X) and J(NU,X), X >= 0, NU >= 0, ACM
  !                 Transactions on Mathematical Software 3, (1977),
  !                 pp. 76-92.
  !               F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !***
  ! **Routines called:**  ASYIK, I1MACH, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4._SP format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH, I1MACH
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alpha, X
  REAL(SP), INTENT(OUT) :: Y(N)
  INTEGER :: i, ialp, in, is, i1, k, kk, km, kt, nn, ns
  REAL(SP) :: ain, ak, akm, ans, ap, arg, atol, tolln, dfn, dtm, dx, earg, elim, &
    etx, flgik, fn, fnf, fni, fnp1, fnu, gln, ra, s, sx, sxo2, s1, s2, t, ta, tb, &
    temp(3), tfn, tm, tol, trx, t2, xo2, xo2l, z
  REAL(SP), PARAMETER :: rttpi = 3.98942280401433E-01_SP
  INTEGER, PARAMETER :: inlim = 80
  !* FIRST EXECUTABLE STATEMENT  BESI
  Nz = 0
  kt = 1
  !     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
  !     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
  ra = R1MACH(3)
  tol = MAX(ra,1.0E-15_SP)
  i1 = -I1MACH(12)
  gln = R1MACH(5)
  elim = 2.303_SP*(i1*gln-3._SP)
  !     TOLLN = -LN(TOL)
  i1 = I1MACH(11) + 1
  tolln = 2.303_SP*gln*i1
  tolln = MIN(tolln,34.5388E0_SP)
  IF( N<1 ) THEN
    ERROR STOP 'BESI : N LESS THAN ONE.'
  ELSEIF( N==1 ) THEN
    kt = 2
  END IF
  nn = N
  IF( Kode<1 .OR. Kode>2 ) THEN
    ERROR STOP 'BESI : SCALING OPTION, KODE, NOT 1 OR 2.'
  ELSEIF( X<0 ) THEN
    ERROR STOP 'BESI : X < 0'
  ELSEIF( Alpha<0._SP ) THEN
    ERROR STOP 'BESI : ORDER, ALPHA < 0 '
  ELSEIF( X==0 ) THEN
    IF( Alpha==0 ) THEN
      Y(1) = 1._SP
      IF( N==1 ) RETURN
      i1 = 2
    ELSE
      i1 = 1
    END IF
    DO i = i1, N
      Y(i) = 0._SP
    END DO
    RETURN
  ELSE
    ialp = INT(Alpha)
    fni = ialp + N - 1
    fnf = Alpha - ialp
    dfn = fni + fnf
    fnu = dfn
    in = 0
    xo2 = X*0.5_SP
    sxo2 = xo2*xo2
    etx = Kode - 1
    sx = etx*X
    !
    !     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
    !     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
    !     APPLIED.
    !
    IF( sxo2<=(fnu+1._SP) ) THEN
      fn = fnu
      fnp1 = fn + 1._SP
      xo2l = LOG(xo2)
      is = kt
      IF( X<=0.5_SP ) GOTO 500
      ns = 0
    ELSEIF( X<=12._SP ) THEN
      xo2l = LOG(xo2)
      ns = INT(sxo2-fnu)
    ELSE
      fn = 0.55_SP*fnu*fnu
      fn = MAX(17._SP,fn)
      IF( X>=fn ) THEN
        !
        !     ASYMPTOTIC EXPANSION FOR X TO INFINITY
        !
        earg = rttpi/SQRT(X)
        IF( Kode==2 ) GOTO 1000
        IF( X>elim ) THEN
          ! CALL XERMSG('BESI : OVERFLOW, X TOO LARGE FOR KODE = 1.',6,1)
          RETURN
        ELSE
          earg = earg*EXP(X)
          GOTO 1000
        END IF
      ELSE
        ans = MAX(36._SP-fnu,0._SP)
        ns = INT(ans)
        fni = fni + ns
        dfn = fni + fnf
        fn = dfn
        is = kt
        km = N - 1 + ns
        IF( km>0 ) is = 3
        !
        !     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
        !
        IF( Kode==2 ) GOTO 100
        IF( Alpha<1._SP ) THEN
          IF( X<=elim ) GOTO 100
          ERROR STOP 'BESI : OVERFLOW, X TOO LARGE FOR KODE = 1.'
        ELSE
          z = X/Alpha
          ra = SQRT(1._SP+z*z)
          gln = LOG((1._SP+ra)/z)
          t = ra*(1._SP-etx) + etx/(z+ra)
          arg = Alpha*(t-gln)
          IF( arg>elim ) THEN
            ! CALL XERMSG('BESI : OVERFLOW, X TOO LARGE FOR KODE = 1.',6,1)
            RETURN
          ELSE
            IF( km/=0 ) GOTO 100
            GOTO 200
          END IF
        END IF
      END IF
    END IF
    fni = fni + ns
    dfn = fni + fnf
    fn = dfn
    fnp1 = fn + 1._SP
    is = kt
    IF( N-1+ns>0 ) is = 3
    GOTO 500
  END IF
  !
  !     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
  !
  100  z = X/fn
  ra = SQRT(1._SP+z*z)
  gln = LOG((1._SP+ra)/z)
  t = ra*(1._SP-etx) + etx/(z+ra)
  arg = fn*(t-gln)
  200 CONTINUE
  IF( arg>=(-elim) ) GOTO 400
  !
  !     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
  !
  Y(nn) = 0._SP
  nn = nn - 1
  fni = fni - 1._SP
  dfn = fni + fnf
  fn = dfn
  IF( nn<1 ) GOTO 800
  IF( nn==1 ) THEN
    kt = 2
    is = 2
  END IF
  GOTO 100
  300  is = 2
  fni = fni - 1._SP
  dfn = fni + fnf
  fn = dfn
  IF( i1==2 ) THEN
    !
    !     BACKWARD RECURSION SECTION
    !
    Nz = N - nn
    GOTO 900
  ELSE
    z = X/fn
    ra = SQRT(1._SP+z*z)
    gln = LOG((1._SP+ra)/z)
    t = ra*(1._SP-etx) + etx/(z+ra)
    arg = fn*(t-gln)
  END IF
  400  i1 = ABS(3-is)
  i1 = MAX(i1,1)
  flgik = 1._SP
  CALL ASYIK(X,fn,Kode,flgik,ra,arg,i1,temp(is))
  SELECT CASE (is)
    CASE (1)
      GOTO 300
    CASE (2)
      Nz = N - nn
      GOTO 900
    CASE (3)
      !     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      t = 1._SP/(fn*ra)
      ain = tolln/(gln+SQRT(gln*gln+t*tolln)) + 1.5_SP
      in = INT(ain)
      IF( in<=inlim ) GOTO 1200
      !
      !     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
      !
      IF( km/=0 ) THEN
        temp(1) = temp(3)
        in = ns
        kt = 1
        i1 = 0
        GOTO 300
      ELSE
        Y(1) = temp(3)
        RETURN
      END IF
  END SELECT
  !
  !     SERIES FOR (X/2)**2<=NU+1
  !
  500  gln = LOG_GAMMA(fnp1)
  arg = fn*xo2l - gln - sx
  IF( arg<(-elim) ) GOTO 700
  earg = EXP(arg)
  600  s = 1._SP
  IF( X>=tol ) THEN
    ak = 3._SP
    t2 = 1._SP
    t = 1._SP
    s1 = fn
    DO k = 1, 17
      s2 = t2 + s1
      t = t*sxo2/s2
      s = s + t
      IF( ABS(t)<tol ) EXIT
      t2 = t2 + ak
      ak = ak + 2._SP
      s1 = s1 + fn
    END DO
  END IF
  temp(is) = s*earg
  SELECT CASE (is)
    CASE (2)
      Nz = N - nn
      GOTO 900
    CASE (3)
      !
      !     BACKWARD RECURSION WITH NORMALIZATION BY
      !     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
      !
      !     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      akm = MAX(3._SP-fn,0._SP)
      km = INT(akm)
      tfn = fn + km
      ta = (gln+tfn-0.9189385332_SP-0.0833333333_SP/tfn)/(tfn+0.5E0_SP)
      ta = xo2l - ta
      tb = -(1._SP-1._SP/tfn)/tfn
      ain = tolln/(-ta+SQRT(ta*ta-tolln*tb)) + 1.5_SP
      in = INT(ain)
      in = in + km
      GOTO 1200
    CASE DEFAULT
      earg = earg*fn/xo2
      fni = fni - 1._SP
      dfn = fni + fnf
      fn = dfn
      is = 2
      GOTO 600
  END SELECT
  700  Y(nn) = 0._SP
  nn = nn - 1
  fnp1 = fn
  fni = fni - 1._SP
  dfn = fni + fnf
  fn = dfn
  IF( nn<1 ) GOTO 800
  IF( nn==1 ) THEN
    kt = 2
    is = 2
  END IF
  IF( sxo2>fnp1 ) GOTO 100
  arg = arg - xo2l + LOG(fnp1)
  IF( arg>=(-elim) ) GOTO 500
  GOTO 700
  800  Nz = N - nn
  RETURN
  900 CONTINUE
  IF( kt==2 ) THEN
    Y(1) = temp(2)
    RETURN
  ELSE
    s1 = temp(1)
    s2 = temp(2)
    trx = 2._SP/X
    dtm = fni
    tm = (dtm+fnf)*trx
    IF( in==0 ) THEN
      !     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(nn) = s1
      Y(nn-1) = s2
      IF( nn==2 ) RETURN
    ELSE
      !     BACKWARD RECUR TO INDEX ALPHA+NN-1
      DO i = 1, in
        s = s2
        s2 = tm*s2 + s1
        s1 = s
        dtm = dtm - 1._SP
        tm = (dtm+fnf)*trx
      END DO
      Y(nn) = s1
      IF( nn==1 ) RETURN
      Y(nn-1) = s2
      IF( nn==2 ) RETURN
    END IF
    k = nn + 1
    DO i = 3, nn
      k = k - 1
      Y(k-2) = tm*Y(k-1) + Y(k)
      dtm = dtm - 1._SP
      tm = (dtm+fnf)*trx
    END DO
    RETURN
  END IF
  1000 etx = 8._SP*X
  is = kt
  in = 0
  fn = fnu
  1100 dx = fni + fni
  tm = 0._SP
  IF( fni/=0._SP .OR. ABS(fnf)>=tol ) tm = 4._SP*fnf*(fni+fni+fnf)
  dtm = dx*dx
  s1 = etx
  trx = dtm - 1._SP
  dx = -(trx+tm)/etx
  t = dx
  s = 1._SP + dx
  atol = tol*ABS(s)
  s2 = 1._SP
  ak = 8._SP
  DO k = 1, 25
    s1 = s1 + etx
    s2 = s2 + ak
    dx = dtm - s2
    ap = dx + tm
    t = -t*ap/s1
    s = s + t
    IF( ABS(t)<=atol ) EXIT
    ak = ak + 8._SP
  END DO
  temp(is) = s*earg
  IF( is==2 ) GOTO 900
  is = 2
  fni = fni - 1._SP
  dfn = fni + fnf
  fn = dfn
  GOTO 1100
  1200 trx = 2._SP/X
  dtm = fni + in
  tm = (dtm+fnf)*trx
  ta = 0._SP
  tb = tol
  kk = 1
  DO
    !
    !     BACKWARD RECUR UNINDEXED
    !
    DO i = 1, in
      s = tb
      tb = tm*tb + ta
      ta = s
      dtm = dtm - 1._SP
      tm = (dtm+fnf)*trx
    END DO
    !     NORMALIZATION
    IF( kk/=1 ) EXIT
    ta = (ta/tb)*temp(3)
    tb = temp(3)
    kk = 2
    in = ns
    IF( ns==0 ) EXIT
  END DO
  Y(nn) = tb
  Nz = N - nn
  IF( nn==1 ) RETURN
  tb = tm*tb + ta
  k = nn - 1
  Y(k) = tb
  IF( nn==2 ) RETURN
  dtm = dtm - 1._SP
  tm = (dtm+fnf)*trx
  km = k - 1
  !
  !     BACKWARD RECUR INDEXED
  !
  DO i = 1, km
    Y(k-1) = tm*Y(k) + Y(k+1)
    dtm = dtm - 1._SP
    tm = (dtm+fnf)*trx
    k = k - 1
  END DO

  RETURN
END SUBROUTINE BESI