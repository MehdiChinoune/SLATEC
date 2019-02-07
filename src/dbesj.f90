!*==DBESJ.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBESJ
SUBROUTINE DBESJ(X,Alpha,N,Y,Nz)
  IMPLICIT NONE
  !*--DBESJ5
  !*** Start of declarations inserted by SPAG
  REAL DJAIRY
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBESJ
  !***PURPOSE  Compute an N member sequence of J Bessel functions
  !            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA
  !            and X.
  !***LIBRARY   SLATEC
  !***CATEGORY  C10A3
  !***TYPE      DOUBLE PRECISION (BESJ-S, DBESJ-D)
  !***KEYWORDS  J BESSEL FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Amos, D. E., (SNLA)
  !           Daniel, S. L., (SNLA)
  !           Weston, M. K., (SNLA)
  !***DESCRIPTION
  !
  !     Abstract  **** a double precision routine ****
  !         DBESJ computes an N member sequence of J Bessel functions
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
  !         when it is efficient to do so. Leading terms of the series and
  !         uniform expansion are tested for underflow.  If a sequence is
  !         requested and the last member would underflow, the result is
  !         set to zero and the next lower order tried, etc., until a
  !         member comes on scale or all members are set to zero.
  !         Overflow cannot occur.
  !
  !         The maximum number of significant digits obtainable
  !         is the smaller of 14 and the number of digits carried in
  !         double precision arithmetic.
  !
  !     Description of Arguments
  !
  !         Input      X,ALPHA are double precision
  !           X      - X .GE. 0.0D0
  !           ALPHA  - order of first member of the sequence,
  !                    ALPHA .GE. 0.0D0
  !           N      - number of members in the sequence, N .GE. 1
  !
  !         Output     Y is double precision
  !           Y      - a vector whose first N components contain
  !                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
  !           NZ     - number of components of Y set to zero due to
  !                    underflow,
  !                    NZ=0   , normal return, computation completed
  !                    NZ .NE. 0, last NZ components of Y set to zero,
  !                             Y(K)=0.0D0, K=N-NZ+1,...,N.
  !
  !     Error Conditions
  !         Improper input arguments - a fatal error
  !         Underflow  - a non-fatal error (NZ .NE. 0)
  !
  !***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
  !                 subroutines IBESS and JBESS for Bessel functions
  !                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
  !                 Transactions on Mathematical Software 3, (1977),
  !                 pp. 76-92.
  !               F. W. J. Olver, Tables of Bessel Functions of Moderate
  !                 or Large Orders, NPL Mathematical Tables 6, Her
  !                 Majesty's Stationery Office, London, 1962.
  !***ROUTINES CALLED  D1MACH, DASYJY, DJAIRY, DLNGAM, I1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   750101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DBESJ
  EXTERNAL DJAIRY
  INTEGER i , ialp , idalp , iflw , in , inlim , is , i1 , i2 , k , kk , &
    km , kt , N , nn , ns , Nz
  INTEGER I1MACH
  DOUBLE PRECISION ak , akm , Alpha , ans , ap , arg , coef , dalpha , dfn , &
    dtm , earg , elim1 , etx , fidal , flgjy , fn , fnf , &
    fni , fnp1 , fnu , fnulim , gln , pdf , pidt , pp , &
    rden , relb , rttp , rtwo , rtx , rzden , s , sa , sb , &
    sxo2 , s1 , s2 , t , ta , tau , tb , temp , tfn , tm , &
    tol , tolln , trx , tx , t1 , t2 , wk , X , xo2 , xo2l , &
    Y , slim , rtol
  SAVE rtwo , pdf , rttp , pidt , pp , inlim , fnulim
  DOUBLE PRECISION D1MACH , DLNGAM
  DIMENSION Y(*) , temp(3) , fnulim(2) , pp(4) , wk(7)
  DATA rtwo , pdf , rttp , pidt/1.34839972492648D+00 , &
    7.85398163397448D-01 , 7.97884560802865D-01 , 1.57079632679490D+00/
  DATA pp(1) , pp(2) , pp(3) , pp(4)/8.72909153935547D+00 , &
    2.65693932265030D-01 , 1.24578576865586D-01 , 7.70133747430388D-04/
  DATA inlim/150/
  DATA fnulim(1) , fnulim(2)/100.0D0 , 60.0D0/
  !***FIRST EXECUTABLE STATEMENT  DBESJ
  Nz = 0
  kt = 1
  ns = 0
  !     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
  !     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
  ta = D1MACH(3)
  tol = MAX(ta,1.0D-15)
  i1 = I1MACH(14) + 1
  i2 = I1MACH(15)
  tb = D1MACH(5)
  elim1 = -2.303D0*(i2*tb+3.0D0)
  rtol = 1.0D0/tol
  slim = D1MACH(1)*rtol*1.0D+3
  !     TOLLN = -LN(TOL)
  tolln = 2.303D0*tb*i1
  tolln = MIN(tolln,34.5388D0)
  IF ( N<1 ) THEN
    CALL XERMSG('SLATEC','DBESJ','N LESS THAN ONE.',2,1)
    RETURN
  ELSEIF ( N==1 ) THEN
    kt = 2
  ENDIF
  nn = N
  IF ( X<0 ) THEN
    CALL XERMSG('SLATEC','DBESJ','X LESS THAN ZERO.',2,1)
    GOTO 99999
  ELSEIF ( X==0 ) THEN
    IF ( Alpha<0 ) GOTO 1200
    IF ( Alpha==0 ) THEN
      Y(1) = 1.0D0
      IF ( N==1 ) RETURN
      i1 = 2
    ELSE
      i1 = 1
    ENDIF
    DO i = i1 , N
      Y(i) = 0.0D0
    ENDDO
    RETURN
  ELSE
    IF ( Alpha<0.0D0 ) GOTO 1200
    !
    ialp = INT(Alpha)
    fni = ialp + N - 1
    fnf = Alpha - ialp
    dfn = fni + fnf
    fnu = dfn
    xo2 = X*0.5D0
    sxo2 = xo2*xo2
    !
    !     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
    !     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
    !     APPLIED.
    !
    IF ( sxo2<=(fnu+1.0D0) ) THEN
      fn = fnu
      fnp1 = fn + 1.0D0
      xo2l = LOG(xo2)
      is = kt
      IF ( X<=0.50D0 ) GOTO 200
      ns = 0
    ELSE
      ta = MAX(20.0D0,fnu)
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
          in = INT(Alpha-tau+2.0D0)
          IF ( in<=0 ) THEN
            idalp = ialp
            in = 0
          ELSE
            idalp = ialp - in - 1
            kt = 1
          ENDIF
          is = kt
          fidal = idalp
          dalpha = fidal + fnf
          arg = X - pidt*dalpha - pdf
          sa = SIN(arg)
          sb = COS(arg)
          coef = rttp/rtx
          etx = 8.0D0*X
          GOTO 800
        ELSE
          fn = fnu
          is = kt
          GOTO 100
        ENDIF
      ELSEIF ( X>12.0D0 ) THEN
        ans = MAX(36.0D0-fnu,0.0D0)
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
      ENDIF
    ENDIF
    fni = fni + ns
    dfn = fni + fnf
    fn = dfn
    fnp1 = fn + 1.0D0
    is = kt
    IF ( N-1+ns>0 ) is = 3
    GOTO 200
  ENDIF
  100  DO
  !
  !     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
  !
  i1 = ABS(3-is)
  i1 = MAX(i1,1)
  flgjy = 1.0D0
  CALL DASYJY(DJAIRY,X,fn,flgjy,i1,temp(is),wk,iflw)
  IF ( iflw/=0 ) THEN
    !
    !     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
    !     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE LARGER
    !     THAN 36. THEREFORE, NS NEE NOT BE TESTED.
    !
    Y(nn) = 0.0D0
    nn = nn - 1
    fni = fni - 1.0D0
    dfn = fni + fnf
    fn = dfn
    IF ( nn<1 ) GOTO 500
    IF ( nn==1 ) THEN
      kt = 2
      is = 2
    ENDIF
  ELSE
    SELECT CASE (is)
      CASE (1)
        EXIT
      CASE (2)
        GOTO 600
      CASE (3)
        !     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
        gln = wk(3) + wk(2)
        IF ( wk(6)>30.0D0 ) THEN
          ta = 0.5D0*tolln/wk(4)
          ta = ((0.0493827160D0*ta-0.1111111111D0)*ta+0.6666666667D0)&
            *ta*wk(6)
          IF ( wk(1)<0.10D0 ) THEN
            tb = (1.259921049D0+(0.1679894730D0+0.0887944358D0*wk(1))*wk(1))&
              /wk(7)
          ELSE
            tb = gln/wk(5)
          ENDIF
        ELSE
          rden = (pp(4)*wk(6)+pp(3))*wk(6) + 1.0D0
          rzden = pp(1) + pp(2)*wk(6)
          ta = rzden/rden
          IF ( wk(1)<0.10D0 ) THEN
            tb = (1.259921049D0+(0.1679894730D0+0.0887944358D0*wk(1))*wk(1))&
              /wk(7)
          ELSE
            tb = gln/wk(5)
          ENDIF
        ENDIF
        in = INT(ta/tb+1.5D0)
        IF ( in<=inlim ) GOTO 900
      CASE DEFAULT
    END SELECT
    temp(1) = temp(3)
    kt = 1
    EXIT
  ENDIF
ENDDO
is = 2
fni = fni - 1.0D0
dfn = fni + fnf
fn = dfn
IF ( i1/=2 ) GOTO 100
GOTO 600
!
!     SERIES FOR (X/2)**2.LE.NU+1
!
200  gln = DLNGAM(fnp1)
arg = fn*xo2l - gln
IF ( arg<(-elim1) ) GOTO 400
earg = EXP(arg)
300  s = 1.0D0
IF ( X>=tol ) THEN
  ak = 3.0D0
  t2 = 1.0D0
  t = 1.0D0
  s1 = fn
  DO k = 1 , 17
    s2 = t2 + s1
    t = -t*sxo2/s2
    s = s + t
    IF ( ABS(t)<tol ) EXIT
    t2 = t2 + ak
    ak = ak + 2.0D0
    s1 = s1 + fn
  ENDDO
ENDIF
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
    akm = MAX(3.0D0-fn,0.0D0)
    km = INT(akm)
    tfn = fn + km
    ta = (gln+tfn-0.9189385332D0-0.0833333333D0/tfn)/(tfn+0.5D0)
    ta = xo2l - ta
    tb = -(1.0D0-1.5D0/tfn)/tfn
    akm = tolln/(-ta+SQRT(ta*ta-tolln*tb)) + 1.5D0
    in = km + INT(akm)
    GOTO 900
  CASE DEFAULT
    earg = earg*fn/xo2
    fni = fni - 1.0D0
    dfn = fni + fnf
    fn = dfn
    is = 2
    GOTO 300
END SELECT
400  Y(nn) = 0.0D0
nn = nn - 1
fnp1 = fn
fni = fni - 1.0D0
dfn = fni + fnf
fn = dfn
IF ( nn<1 ) GOTO 500
IF ( nn==1 ) THEN
  kt = 2
  is = 2
ENDIF
IF ( sxo2>fnp1 ) GOTO 100
arg = arg - xo2l + LOG(fnp1)
IF ( arg>=(-elim1) ) GOTO 200
GOTO 400
500  Nz = N - nn
RETURN
!
!     BACKWARD RECURSION SECTION
!
600  IF ( ns==0 ) THEN
Nz = N - nn
IF ( kt==2 ) GOTO 700
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
Y(nn) = temp(1)
Y(nn-1) = temp(2)
IF ( nn==2 ) RETURN
ENDIF
trx = 2.0D0/X
dtm = fni
tm = (dtm+fnf)*trx
ak = 1.0D0
ta = temp(1)
tb = temp(2)
IF ( ABS(ta)<=slim ) THEN
ta = ta*rtol
tb = tb*rtol
ak = tol
ENDIF
kk = 2
in = ns - 1
IF ( in==0 ) GOTO 1100
IF ( ns/=0 ) GOTO 1000
k = nn - 2
DO i = 3 , nn
s = tb
tb = tm*tb - ta
ta = s
Y(k) = tb*ak
dtm = dtm - 1.0D0
tm = (dtm+fnf)*trx
k = k - 1
ENDDO
RETURN
700  Y(1) = temp(2)
RETURN
800  dtm = fidal + fidal
dtm = dtm*dtm
tm = 0.0D0
IF ( fidal/=0.0D0.OR.ABS(fnf)>=tol ) tm = 4.0D0*fnf*(fidal+fidal+fnf)
trx = dtm - 1.0D0
t2 = (trx+tm)/etx
s2 = t2
relb = tol*ABS(t2)
t1 = etx
s1 = 1.0D0
fn = 1.0D0
ak = 8.0D0
DO k = 1 , 13
t1 = t1 + etx
fn = fn + ak
trx = dtm - fn
ap = trx + tm
t2 = -t2*ap/t1
s1 = s1 + t2
t1 = t1 + etx
ak = ak + 8.0D0
fn = fn + ak
trx = dtm - fn
ap = trx + tm
t2 = t2*ap/t1
s2 = s2 + t2
IF ( ABS(t2)<=relb ) EXIT
ak = ak + 8.0D0
ENDDO
temp(is) = coef*(s1*sb-s2*sa)
IF ( is==2 ) THEN
!
!     FORWARD RECURSION SECTION
!
IF ( kt==2 ) GOTO 700
s1 = temp(1)
s2 = temp(2)
tx = 2.0D0/X
tm = dalpha*tx
IF ( in/=0 ) THEN
  !
  !     FORWARD RECUR TO INDEX ALPHA
  !
  DO i = 1 , in
    s = s2
    s2 = tm*s2 - s1
    tm = tm + tx
    s1 = s
  ENDDO
  IF ( nn==1 ) THEN
    Y(1) = s2
    RETURN
  ELSE
    s = s2
    s2 = tm*s2 - s1
    tm = tm + tx
    s1 = s
  ENDIF
ENDIF
!
!     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
!
Y(1) = s1
Y(2) = s2
IF ( nn==2 ) RETURN
DO i = 3 , nn
  Y(i) = tm*Y(i-1) - Y(i-2)
  tm = tm + tx
ENDDO
RETURN
ELSE
fidal = fidal + 1.0D0
dalpha = fidal + fnf
is = 2
tb = sa
sa = -sb
sb = tb
GOTO 800
ENDIF
900  dtm = fni + in
trx = 2.0D0/X
tm = (dtm+fnf)*trx
ta = 0.0D0
tb = tol
kk = 1
ak = 1.0D0
1000 DO
!
!     BACKWARD RECUR UNINDEXED
!
DO i = 1 , in
s = tb
tb = tm*tb - ta
ta = s
dtm = dtm - 1.0D0
tm = (dtm+fnf)*trx
ENDDO
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
ENDIF
ta = ta*sa
kk = 2
in = ns
IF ( ns==0 ) EXIT
ENDDO
1100 Y(nn) = tb*ak
Nz = N - nn
IF ( nn==1 ) RETURN
k = nn - 1
s = tb
tb = tm*tb - ta
ta = s
Y(k) = tb*ak
IF ( nn==2 ) RETURN
dtm = dtm - 1.0D0
tm = (dtm+fnf)*trx
k = nn - 2
!
!     BACKWARD RECUR INDEXED
!
DO i = 3 , nn
s = tb
tb = tm*tb - ta
ta = s
Y(k) = tb*ak
dtm = dtm - 1.0D0
tm = (dtm+fnf)*trx
k = k - 1
ENDDO
RETURN
!
!
!
1200 CALL XERMSG('SLATEC','DBESJ','ORDER, ALPHA, LESS THAN ZERO.',2,1)
RETURN
99999 END SUBROUTINE DBESJ
