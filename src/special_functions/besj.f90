!** BESJ
SUBROUTINE BESJ(X,Alpha,N,Y,Nz)
  IMPLICIT NONE
  !>
  !***
  !  Compute an N member sequence of J Bessel functions
  !            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA
  !            and X.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C10A3
  !***
  ! **Type:**      SINGLE PRECISION (BESJ-S, DBESJ-D)
  !***
  ! **Keywords:**  J BESSEL FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !           Daniel, S. L., (SNLA)
  !           Weston, M. K., (SNLA)
  !***
  ! **Description:**
  !
  !     Abstract
  !         BESJ computes an N member sequence of J Bessel functions
  !         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
  !         A combination of the power series, the asymptotic expansion
  !         for X to infinity and the uniform asymptotic expansion for
  !         NU to infinity are applied over subdivisions of the (NU,X)
  !         plane.  For values of (NU,X) not covered by one of these
  !         formulae, the order is incremented or decremented by integer
  !         values into a region where one of the formulae apply. Backward
  !         recursion is applied to reduce orders by integer values except
  !         where the entire sequence lies in the oscillatory region.  In
  !         this case forward recursion is stable and values from the
  !         asymptotic expansion for X to infinity start the recursion
  !         when it is efficient to do so.  Leading terms of the series
  !         and uniform expansion are tested for underflow.  If a sequence
  !         is requested and the last member would underflow, the result
  !         is set to zero and the next lower order tried, etc., until a
  !         member comes on scale or all members are set to zero.
  !         Overflow cannot occur.
  !
  !     Description of Arguments
  !
  !         Input
  !           X      - X .GE. 0.0E0
  !           ALPHA  - order of first member of the sequence,
  !                    ALPHA .GE. 0.0E0
  !           N      - number of members in the sequence, N .GE. 1
  !
  !         Output
  !           Y      - a vector whose first  N components contain
  !                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
  !           NZ     - number of components of Y set to zero due to
  !                    underflow,
  !                    NZ=0  , normal return, computation completed
  !                    NZ .NE. 0, last NZ components of Y set to zero,
  !                             Y(K)=0.0E0, K=N-NZ+1,...,N.
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Underflow  - a non-fatal error (NZ .NE. 0)
  !
  !***
  ! **References:**  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
  !                 subroutines IBESS and JBESS for Bessel functions
  !                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
  !                 Transactions on Mathematical Software 3, (1977),
  !                 pp. 76-92.
  !               F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !***
  ! **Routines called:**  ASYJY, I1MACH, JAIRY, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER i, ialp, idalp, iflw, in, is, i1, i2, k, kk, km, kt, N, nn, ns, Nz
  REAL ak, akm, Alpha, ans, ap, arg, coef, dalpha, dfn, dtm, earg, elim1, etx, &
    fidal, flgjy, fn, fnf, fni, fnp1, fnu, gln, rden, relb, rtx, rzden, s, sa, &
    sb, sxo2, s1, s2, t, ta, tau, tb, temp(3), tfn, tm, tol, tolln, trx, tx, t1, &
    t2, wk(7), X, xo2, xo2l, Y(*), rtol, slim
  INTEGER, EXTERNAL :: I1MACH
  REAL, EXTERNAL :: R1MACH
  EXTERNAL :: JAIRY
  REAL, PARAMETER :: rtwo = 1.34839972492648E+00, pdf = 7.85398163397448E-01, &
    rttp = 7.97884560802865E-01, pidt = 1.57079632679490E+00
  REAL, PARAMETER :: pp(4) = [ 8.72909153935547E+00, 2.65693932265030E-01, &
    1.24578576865586E-01, 7.70133747430388E-04 ]
  INTEGER, PARAMETER :: inlim = 150
  REAL, PARAMETER :: fnulim(2) = [ 100.0E0, 60.0E0 ]
  !* FIRST EXECUTABLE STATEMENT  BESJ
  Nz = 0
  kt = 1
  ns = 0
  !     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
  !     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
  ta = R1MACH(3)
  tol = MAX(ta,1.0E-15)
  i1 = I1MACH(11) + 1
  i2 = I1MACH(12)
  tb = R1MACH(5)
  elim1 = -2.303E0*(i2*tb+3.0E0)
  rtol = 1.0E0/tol
  slim = R1MACH(1)*1.0E+3*rtol
  !     TOLLN = -LN(TOL)
  tolln = 2.303E0*tb*i1
  tolln = MIN(tolln,34.5388E0)
  IF ( N<1 ) THEN
    CALL XERMSG('SLATEC','BESJ','N LESS THAN ONE.',2,1)
    RETURN
  ELSEIF ( N==1 ) THEN
    kt = 2
  END IF
  nn = N
  IF ( X<0 ) THEN
    CALL XERMSG('SLATEC','BESJ','X LESS THAN ZERO.',2,1)
    RETURN
  ELSEIF ( X==0 ) THEN
    IF ( Alpha<0 ) GOTO 1200
    IF ( Alpha==0 ) THEN
      Y(1) = 1.0E0
      IF ( N==1 ) RETURN
      i1 = 2
    ELSE
      i1 = 1
    END IF
    DO i = i1, N
      Y(i) = 0.0E0
    END DO
    RETURN
  ELSE
    IF ( Alpha<0.0E0 ) GOTO 1200
    !
    ialp = INT(Alpha)
    fni = ialp + N - 1
    fnf = Alpha - ialp
    dfn = fni + fnf
    fnu = dfn
    xo2 = X*0.5E0
    sxo2 = xo2*xo2
    !
    !     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
    !     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
    !     APPLIED.
    !
    IF ( sxo2<=(fnu+1.0E0) ) THEN
      fn = fnu
      fnp1 = fn + 1.0E0
      xo2l = LOG(xo2)
      is = kt
      IF ( X<=0.50E0 ) GOTO 200
      ns = 0
    ELSE
      ta = MAX(20.0E0,fnu)
      IF ( X>ta ) THEN
        rtx = SQRT(X)
        tau = rtwo*rtx
        ta = tau + fnulim(kt)
        IF ( fnu<=ta ) THEN
          !
          !     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
          !     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
          !     OF THE SEQUENCE IS ALSO IN THE REGION.
          !
          in = INT(Alpha-tau+2.0E0)
          IF ( in<=0 ) THEN
            idalp = ialp
            in = 0
          ELSE
            idalp = ialp - in - 1
            kt = 1
          END IF
          is = kt
          fidal = idalp
          dalpha = fidal + fnf
          arg = X - pidt*dalpha - pdf
          sa = SIN(arg)
          sb = COS(arg)
          coef = rttp/rtx
          etx = 8.0E0*X
          GOTO 800
        ELSE
          fn = fnu
          is = kt
          GOTO 100
        END IF
      ELSEIF ( X>12.0E0 ) THEN
        ans = MAX(36.0E0-fnu,0.0E0)
        ns = INT(ans)
        fni = fni + ns
        dfn = fni + fnf
        fn = dfn
        is = kt
        IF ( N-1+ns>0 ) is = 3
        GOTO 100
      ELSE
        xo2l = LOG(xo2)
        ns = INT(sxo2-fnu) + 1
      END IF
    END IF
    fni = fni + ns
    dfn = fni + fnf
    fn = dfn
    fnp1 = fn + 1.0E0
    is = kt
    IF ( N-1+ns>0 ) is = 3
    GOTO 200
  END IF
  100 CONTINUE
  DO
    !
    !     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
    !
    i1 = ABS(3-is)
    i1 = MAX(i1,1)
    flgjy = 1.0E0
    CALL ASYJY(JAIRY,X,fn,flgjy,i1,temp(is),wk,iflw)
    IF ( iflw/=0 ) THEN
      !
      !     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
      !     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE
      !     LARGER THAN 36. THEREFORE, NS NEED NOT BE CONSIDERED.
      !
      Y(nn) = 0.0E0
      nn = nn - 1
      fni = fni - 1.0E0
      dfn = fni + fnf
      fn = dfn
      IF ( nn<1 ) GOTO 500
      IF ( nn==1 ) THEN
        kt = 2
        is = 2
      END IF
    ELSE
      SELECT CASE (is)
        CASE (1)
          EXIT
        CASE (2)
          GOTO 600
        CASE (3)
          !     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
          gln = wk(3) + wk(2)
          IF ( wk(6)>30.0E0 ) THEN
            ta = 0.5E0*tolln/wk(4)
            ta = ((0.0493827160E0*ta-0.1111111111E0)*ta+0.6666666667E0)*ta*wk(6)
            IF ( wk(1)<0.10E0 ) THEN
              tb = (1.259921049E0+(0.1679894730E0+0.0887944358E0*wk(1))*wk(1))&
                /wk(7)
            ELSE
              tb = gln/wk(5)
            END IF
          ELSE
            rden = (pp(4)*wk(6)+pp(3))*wk(6) + 1.0E0
            rzden = pp(1) + pp(2)*wk(6)
            ta = rzden/rden
            IF ( wk(1)<0.10E0 ) THEN
              tb = (1.259921049E0+(0.1679894730E0+0.0887944358E0*wk(1))*wk(1))&
                /wk(7)
            ELSE
              tb = gln/wk(5)
            END IF
          END IF
          in = INT(ta/tb+1.5E0)
          IF ( in<=inlim ) GOTO 900
        CASE DEFAULT
      END SELECT
      temp(1) = temp(3)
      kt = 1
      EXIT
    END IF
  END DO
  is = 2
  fni = fni - 1.0E0
  dfn = fni + fnf
  fn = dfn
  IF ( i1/=2 ) GOTO 100
  GOTO 600
  !
  !     SERIES FOR (X/2)**2.LE.NU+1
  !
  200  gln = LOG_GAMMA(fnp1)
  arg = fn*xo2l - gln
  IF ( arg<(-elim1) ) GOTO 400
  earg = EXP(arg)
  300  s = 1.0E0
  IF ( X>=tol ) THEN
    ak = 3.0E0
    t2 = 1.0E0
    t = 1.0E0
    s1 = fn
    DO k = 1, 17
      s2 = t2 + s1
      t = -t*sxo2/s2
      s = s + t
      IF ( ABS(t)<tol ) EXIT
      t2 = t2 + ak
      ak = ak + 2.0E0
      s1 = s1 + fn
    END DO
  END IF
  temp(is) = s*earg
  SELECT CASE (is)
    CASE (2)
      GOTO 600
    CASE (3)
      !
      !     BACKWARD RECURSION WITH NORMALIZATION BY
      !     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
      !
      !     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      akm = MAX(3.0E0-fn,0.0E0)
      km = INT(akm)
      tfn = fn + km
      ta = (gln+tfn-0.9189385332E0-0.0833333333E0/tfn)/(tfn+0.5E0)
      ta = xo2l - ta
      tb = -(1.0E0-1.5E0/tfn)/tfn
      akm = tolln/(-ta+SQRT(ta*ta-tolln*tb)) + 1.5E0
      in = km + INT(akm)
      GOTO 900
    CASE DEFAULT
      earg = earg*fn/xo2
      fni = fni - 1.0E0
      dfn = fni + fnf
      fn = dfn
      is = 2
      GOTO 300
  END SELECT
  400  Y(nn) = 0.0E0
  nn = nn - 1
  fnp1 = fn
  fni = fni - 1.0E0
  dfn = fni + fnf
  fn = dfn
  IF ( nn<1 ) GOTO 500
  IF ( nn==1 ) THEN
    kt = 2
    is = 2
  END IF
  IF ( sxo2>fnp1 ) GOTO 100
  arg = arg - xo2l + LOG(fnp1)
  IF ( arg>=(-elim1) ) GOTO 200
  GOTO 400
  500  Nz = N - nn
  RETURN
  !
  !     BACKWARD RECURSION SECTION
  !
  600 CONTINUE
  IF ( ns==0 ) THEN
    Nz = N - nn
    IF ( kt==2 ) GOTO 700
    !     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
    Y(nn) = temp(1)
    Y(nn-1) = temp(2)
    IF ( nn==2 ) RETURN
  END IF
  trx = 2.0E0/X
  dtm = fni
  tm = (dtm+fnf)*trx
  ak = 1.0E0
  ta = temp(1)
  tb = temp(2)
  IF ( ABS(ta)<=slim ) THEN
    ta = ta*rtol
    tb = tb*rtol
    ak = tol
  END IF
  kk = 2
  in = ns - 1
  IF ( in==0 ) GOTO 1100
  IF ( ns/=0 ) GOTO 1000
  k = nn - 2
  DO i = 3, nn
    s = tb
    tb = tm*tb - ta
    ta = s
    Y(k) = tb*ak
    k = k - 1
    dtm = dtm - 1.0E0
    tm = (dtm+fnf)*trx
  END DO
  RETURN
  700  Y(1) = temp(2)
  RETURN
  800  dtm = fidal + fidal
  dtm = dtm*dtm
  tm = 0.0E0
  IF ( fidal/=0.0E0.OR.ABS(fnf)>=tol ) tm = 4.0E0*fnf*(fidal+fidal+fnf)
  trx = dtm - 1.0E0
  t2 = (trx+tm)/etx
  s2 = t2
  relb = tol*ABS(t2)
  t1 = etx
  s1 = 1.0E0
  fn = 1.0E0
  ak = 8.0E0
  DO k = 1, 13
    t1 = t1 + etx
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = -t2*ap/t1
    s1 = s1 + t2
    t1 = t1 + etx
    ak = ak + 8.0E0
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = t2*ap/t1
    s2 = s2 + t2
    IF ( ABS(t2)<=relb ) EXIT
    ak = ak + 8.0E0
  END DO
  temp(is) = coef*(s1*sb-s2*sa)
  IF ( is==2 ) THEN
    !
    !     FORWARD RECURSION SECTION
    !
    IF ( kt==2 ) GOTO 700
    s1 = temp(1)
    s2 = temp(2)
    tx = 2.0E0/X
    tm = dalpha*tx
    IF ( in/=0 ) THEN
      !
      !     FORWARD RECUR TO INDEX ALPHA
      !
      DO i = 1, in
        s = s2
        s2 = tm*s2 - s1
        tm = tm + tx
        s1 = s
      END DO
      IF ( nn==1 ) THEN
        Y(1) = s2
        RETURN
      ELSE
        s = s2
        s2 = tm*s2 - s1
        tm = tm + tx
        s1 = s
      END IF
    END IF
    !
    !     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
    !
    Y(1) = s1
    Y(2) = s2
    IF ( nn==2 ) RETURN
    DO i = 3, nn
      Y(i) = tm*Y(i-1) - Y(i-2)
      tm = tm + tx
    END DO
    RETURN
  ELSE
    fidal = fidal + 1.0E0
    dalpha = fidal + fnf
    is = 2
    tb = sa
    sa = -sb
    sb = tb
    GOTO 800
  END IF
  900  dtm = fni + in
  trx = 2.0E0/X
  tm = (dtm+fnf)*trx
  ta = 0.0E0
  tb = tol
  kk = 1
  ak = 1.0E0
  1000 CONTINUE
  DO
    !
    !     BACKWARD RECUR UNINDEXED AND SCALE WHEN MAGNITUDES ARE CLOSE TO
    !     UNDERFLOW LIMITS (LESS THAN SLIM=R1MACH(1)*1.0E+3/TOL)
    !
    DO i = 1, in
      s = tb
      tb = tm*tb - ta
      ta = s
      dtm = dtm - 1.0E0
      tm = (dtm+fnf)*trx
    END DO
    !     NORMALIZATION
    IF ( kk/=1 ) EXIT
    s = temp(3)
    sa = ta/tb
    ta = s
    tb = s
    IF ( ABS(s)<=slim ) THEN
      ta = ta*rtol
      tb = tb*rtol
      ak = tol
    END IF
    ta = ta*sa
    kk = 2
    in = ns
    IF ( ns==0 ) EXIT
  END DO
  1100 Y(nn) = tb*ak
  Nz = N - nn
  IF ( nn==1 ) RETURN
  k = nn - 1
  s = tb
  tb = tm*tb - ta
  ta = s
  Y(k) = tb*ak
  IF ( nn==2 ) RETURN
  dtm = dtm - 1.0E0
  tm = (dtm+fnf)*trx
  k = nn - 2
  !
  !     BACKWARD RECUR INDEXED
  !
  DO i = 3, nn
    s = tb
    tb = tm*tb - ta
    ta = s
    Y(k) = tb*ak
    dtm = dtm - 1.0E0
    tm = (dtm+fnf)*trx
    k = k - 1
  END DO
  RETURN
  !
  !
  !
  1200 CALL XERMSG('SLATEC','BESJ','ORDER, ALPHA, LESS THAN ZERO.',2,1)
  RETURN
END SUBROUTINE BESJ
