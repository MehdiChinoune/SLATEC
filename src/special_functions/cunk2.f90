!DECK CUNK2
SUBROUTINE CUNK2(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CUNK2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUNK2-A, ZUNK2-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
  !     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
  !     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
  !     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
  !     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
  !     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
  !     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
  !     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
  !
  !***SEE ALSO  CBESK
  !***ROUTINES CALLED  CAIRY, CS1S2, CUCHK, CUNHJ, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CUNK2
  COMPLEX ai, arg, asum, bsum, cfn, ci, cip, ck, cone, crsc, cr1, &
    cr2, cs, cscl, csgn, cspn, csr, css, cy, czero, c1, c2, &
    dai, phi, rz, s1, s2, Y, Z, zb, zeta1, zeta2, zn, zr, &
    phid, argd, zeta1d, zeta2d, asumd, bsumd
  REAL aarg, aic, Alim, ang, aphi, asc, ascle, bry, car, cpn, &
    c2i, c2m, c2r, Elim, fmr, fn, fnf, Fnu, hpi, pi, rs1, &
    sar, sgn, spn, Tol, x, yy, R1MACH
  INTEGER i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, &
    kk, Kode, Mr, N, nai, ndai, nw, Nz, idum, j, ipard, ic
  DIMENSION bry(3), Y(N), asum(2), bsum(2), phi(2), arg(2), zeta1(2), &
    zeta2(2), cy(2), cip(4), css(3), csr(3)
  DATA czero, cone, ci, cr1, cr2/(0.0E0,0.0E0), (1.0E0,0.0E0), &
    (0.0E0,1.0E0), (1.0E0,1.73205080756887729E0), &
    (-0.5E0,-8.66025403784438647E-01)/
  DATA hpi, pi, aic/1.57079632679489662E+00, 3.14159265358979324E+00, &
    1.26551212348464539E+00/
  DATA cip(1), cip(2), cip(3), cip(4)/(1.0E0,0.0E0), (0.0E0,-1.0E0), &
    (-1.0E0,0.0E0), (0.0E0,1.0E0)/
  !***FIRST EXECUTABLE STATEMENT  CUNK2
  kdflg = 1
  Nz = 0
  !-----------------------------------------------------------------------
  !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
  !     THE UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  cscl = CMPLX(1.0E0/Tol,0.0E0)
  crsc = CMPLX(Tol,0.0E0)
  css(1) = cscl
  css(2) = cone
  css(3) = crsc
  csr(1) = crsc
  csr(2) = cone
  csr(3) = cscl
  bry(1) = 1.0E+3*R1MACH(1)/Tol
  bry(2) = 1.0E0/bry(1)
  bry(3) = R1MACH(2)
  x = REAL(Z)
  zr = Z
  IF ( x<0.0E0 ) zr = -Z
  yy = AIMAG(zr)
  zn = -zr*ci
  zb = zr
  inu = INT( Fnu )
  fnf = Fnu - inu
  ang = -hpi*fnf
  car = COS(ang)
  sar = SIN(ang)
  cpn = -hpi*car
  spn = -hpi*sar
  c2 = CMPLX(-spn,cpn)
  kk = MOD(inu,4) + 1
  cs = cr1*c2*cip(kk)
  IF ( yy<=0.0E0 ) THEN
    zn = CONJG(-zn)
    zb = CONJG(zb)
  ENDIF
  !-----------------------------------------------------------------------
  !     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
  !     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
  !     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
  !-----------------------------------------------------------------------
  j = 2
  DO i = 1, N
    !-----------------------------------------------------------------------
    !     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
    !-----------------------------------------------------------------------
    j = 3 - j
    fn = Fnu + (i-1)
    CALL CUNHJ(zn,fn,0,Tol,phi(j),arg(j),zeta1(j),zeta2(j),asum(j),bsum(j))
    IF ( Kode==1 ) THEN
      s1 = zeta1(j) - zeta2(j)
    ELSE
      cfn = CMPLX(fn,0.0E0)
      s1 = zeta1(j) - cfn*(cfn/(zb+zeta2(j)))
    ENDIF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1)
    IF ( ABS(rs1)<=Elim ) THEN
      IF ( kdflg==1 ) kflag = 2
      IF ( ABS(rs1)>=Alim ) THEN
        !-----------------------------------------------------------------------
        !     REFINE  TEST AND SCALE
        !-----------------------------------------------------------------------
        aphi = ABS(phi(j))
        aarg = ABS(arg(j))
        rs1 = rs1 + ALOG(aphi) - 0.25E0*ALOG(aarg) - aic
        IF ( ABS(rs1)>Elim ) GOTO 50
        IF ( kdflg==1 ) kflag = 1
        IF ( rs1>=0.0E0 ) THEN
          IF ( kdflg==1 ) kflag = 3
        ENDIF
      ENDIF
      !-----------------------------------------------------------------------
      !     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
      !     EXPONENT EXTREMES
      !-----------------------------------------------------------------------
      c2 = arg(j)*cr2
      CALL CAIRY(c2,0,2,ai,nai,idum)
      CALL CAIRY(c2,1,2,dai,ndai,idum)
      s2 = cs*phi(j)*(ai*asum(j)+cr2*dai*bsum(j))
      c2r = REAL(s1)
      c2i = AIMAG(s1)
      c2m = EXP(c2r)*REAL(css(kflag))
      s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
      s2 = s2*s1
      IF ( kflag==1 ) THEN
        CALL CUCHK(s2,nw,bry(1),Tol)
        IF ( nw/=0 ) GOTO 50
      ENDIF
      IF ( yy<=0.0E0 ) s2 = CONJG(s2)
      cy(kdflg) = s2
      Y(i) = s2*csr(kflag)
      cs = -ci*cs
      IF ( kdflg==2 ) GOTO 100
      kdflg = 2
      CYCLE
    ENDIF
    50     IF ( rs1>0.0E0 ) GOTO 600
    !-----------------------------------------------------------------------
    !     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
    !-----------------------------------------------------------------------
    IF ( x<0.0E0 ) GOTO 600
    kdflg = 1
    Y(i) = czero
    cs = -ci*cs
    Nz = Nz + 1
    IF ( i/=1 ) THEN
      IF ( Y(i-1)/=czero ) THEN
        Y(i-1) = czero
        Nz = Nz + 1
      ENDIF
    ENDIF
  ENDDO
  i = N
  100  rz = CMPLX(2.0E0,0.0E0)/zr
  ck = CMPLX(fn,0.0E0)*rz
  ib = i + 1
  IF ( N<ib ) GOTO 300
  !-----------------------------------------------------------------------
  !     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
  !     ON UNDERFLOW
  !-----------------------------------------------------------------------
  fn = Fnu + (N-1)
  ipard = 1
  IF ( Mr/=0 ) ipard = 0
  CALL CUNHJ(zn,fn,ipard,Tol,phid,argd,zeta1d,zeta2d,asumd,bsumd)
  IF ( Kode==1 ) THEN
    s1 = zeta1d - zeta2d
  ELSE
    cfn = CMPLX(fn,0.0E0)
    s1 = zeta1d - cfn*(cfn/(zb+zeta2d))
  ENDIF
  rs1 = REAL(s1)
  IF ( ABS(rs1)<=Elim ) THEN
    IF ( ABS(rs1)<Alim ) GOTO 200
    !-----------------------------------------------------------------------
    !     REFINE ESTIMATE AND TEST
    !-----------------------------------------------------------------------
    aphi = ABS(phid)
    aarg = ABS(argd)
    rs1 = rs1 + ALOG(aphi) - 0.25E0*ALOG(aarg) - aic
    IF ( ABS(rs1)<Elim ) GOTO 200
  ENDIF
  IF ( rs1>0.0E0 ) GOTO 600
  !-----------------------------------------------------------------------
  !     FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
  !-----------------------------------------------------------------------
  IF ( x<0.0E0 ) GOTO 600
  Nz = N
  DO i = 1, N
    Y(i) = czero
  ENDDO
  RETURN
  !-----------------------------------------------------------------------
  !     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
  !-----------------------------------------------------------------------
  200  s1 = cy(1)
  s2 = cy(2)
  c1 = csr(kflag)
  ascle = bry(kflag)
  DO i = ib, N
    c2 = s2
    s2 = ck*s2 + s1
    s1 = c2
    ck = ck + rz
    c2 = s2*c1
    Y(i) = c2
    IF ( kflag<3 ) THEN
      c2r = REAL(c2)
      c2i = AIMAG(c2)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF ( c2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*c1
        s2 = c2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        c1 = csr(kflag)
      ENDIF
    ENDIF
  ENDDO
  300 CONTINUE
  IF ( Mr==0 ) RETURN
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
  !-----------------------------------------------------------------------
  Nz = 0
  fmr = Mr
  sgn = -SIGN(pi,fmr)
  !-----------------------------------------------------------------------
  !     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
  !-----------------------------------------------------------------------
  csgn = CMPLX(0.0E0,sgn)
  IF ( yy<=0.0E0 ) csgn = CONJG(csgn)
  ifn = inu + N - 1
  ang = fnf*sgn
  cpn = COS(ang)
  spn = SIN(ang)
  cspn = CMPLX(cpn,spn)
  IF ( MOD(ifn,2)==1 ) cspn = -cspn
  !-----------------------------------------------------------------------
  !     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
  !     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
  !     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
  !     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
  !-----------------------------------------------------------------------
  cs = CMPLX(car,-sar)*csgn
  in = MOD(ifn,4) + 1
  c2 = cip(in)
  cs = cs*CONJG(c2)
  asc = bry(1)
  kk = N
  kdflg = 1
  ib = ib - 1
  ic = ib - 1
  iuf = 0
  DO k = 1, N
    !-----------------------------------------------------------------------
    !     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
    !     FUNCTION ABOVE
    !-----------------------------------------------------------------------
    fn = Fnu + (kk-1)
    IF ( N>2 ) THEN
      IF ( (kk==N).AND.(ib<N) ) GOTO 350
      IF ( (kk/=ib).AND.(kk/=ic) ) THEN
        CALL CUNHJ(zn,fn,0,Tol,phid,argd,zeta1d,zeta2d,asumd,bsumd)
        GOTO 350
      ENDIF
    ENDIF
    phid = phi(j)
    argd = arg(j)
    zeta1d = zeta1(j)
    zeta2d = zeta2(j)
    asumd = asum(j)
    bsumd = bsum(j)
    j = 3 - j
    350    IF ( Kode==1 ) THEN
    s1 = -zeta1d + zeta2d
  ELSE
    cfn = CMPLX(fn,0.0E0)
    s1 = -zeta1d + cfn*(cfn/(zb+zeta2d))
  ENDIF
  !-----------------------------------------------------------------------
  !     TEST FOR UNDERFLOW AND OVERFLOW
  !-----------------------------------------------------------------------
  rs1 = REAL(s1)
  IF ( ABS(rs1)>Elim ) GOTO 450
  IF ( kdflg==1 ) iflag = 2
  IF ( ABS(rs1)>=Alim ) THEN
    !-----------------------------------------------------------------------
    !     REFINE  TEST AND SCALE
    !-----------------------------------------------------------------------
    aphi = ABS(phid)
    aarg = ABS(argd)
    rs1 = rs1 + ALOG(aphi) - 0.25E0*ALOG(aarg) - aic
    IF ( ABS(rs1)>Elim ) GOTO 450
    IF ( kdflg==1 ) iflag = 1
    IF ( rs1>=0.0E0 ) THEN
      IF ( kdflg==1 ) iflag = 3
    ENDIF
  ENDIF
  CALL CAIRY(argd,0,2,ai,nai,idum)
  CALL CAIRY(argd,1,2,dai,ndai,idum)
  s2 = cs*phid*(ai*asumd+dai*bsumd)
  c2r = REAL(s1)
  c2i = AIMAG(s1)
  c2m = EXP(c2r)*REAL(css(iflag))
  s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
  s2 = s2*s1
  IF ( iflag==1 ) THEN
    CALL CUCHK(s2,nw,bry(1),Tol)
    IF ( nw/=0 ) s2 = CMPLX(0.0E0,0.0E0)
  ENDIF
  400    IF ( yy<=0.0E0 ) s2 = CONJG(s2)
  cy(kdflg) = s2
  c2 = s2
  s2 = s2*csr(iflag)
  !-----------------------------------------------------------------------
  !     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
  !-----------------------------------------------------------------------
  s1 = Y(kk)
  IF ( Kode/=1 ) THEN
    CALL CS1S2(zr,s1,s2,nw,asc,Alim,iuf)
    Nz = Nz + nw
  ENDIF
  Y(kk) = s1*cspn + s2
  kk = kk - 1
  cspn = -cspn
  cs = -cs*ci
  IF ( c2/=czero ) THEN
    IF ( kdflg==2 ) GOTO 500
    kdflg = 2
    CYCLE
  ELSE
    kdflg = 1
    CYCLE
  ENDIF
  450    IF ( rs1>0.0E0 ) GOTO 600
  s2 = czero
  GOTO 400
  ENDDO
  k = N
  500  il = N - k
  IF ( il==0 ) RETURN
  !-----------------------------------------------------------------------
  !     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
  !     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
  !     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
  !-----------------------------------------------------------------------
  s1 = cy(1)
  s2 = cy(2)
  cs = csr(iflag)
  ascle = bry(iflag)
  fn = inu + il
  DO i = 1, il
    c2 = s2
    s2 = s1 + CMPLX(fn+fnf,0.0E0)*rz*s2
    s1 = c2
    fn = fn - 1.0E0
    c2 = s2*cs
    ck = c2
    c1 = Y(kk)
    IF ( Kode/=1 ) THEN
      CALL CS1S2(zr,c1,c2,nw,asc,Alim,iuf)
      Nz = Nz + nw
    ENDIF
    Y(kk) = c1*cspn + c2
    kk = kk - 1
    cspn = -cspn
    IF ( iflag<3 ) THEN
      c2r = REAL(ck)
      c2i = AIMAG(ck)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF ( c2m>ascle ) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1*cs
        s2 = ck
        s1 = s1*css(iflag)
        s2 = s2*css(iflag)
        cs = csr(iflag)
      ENDIF
    ENDIF
  ENDDO
  RETURN
  600  Nz = -1
END SUBROUTINE CUNK2
