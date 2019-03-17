!DECK ZUNK2
SUBROUTINE ZUNK2(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ZUNK2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUNK2-A, ZUNK2-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
  !     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
  !     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
  !     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
  !     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
  !     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
  !     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
  !     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
  !
  !***SEE ALSO  ZBESK
  !***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZS1S2, ZUCHK, ZUNHJ
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZUNK2
  !     COMPLEX AI,ARG,ARGD,ASUM,ASUMD,BSUM,BSUMD,CFN,CI,CIP,CK,CONE,CRSC,
  !    *CR1,CR2,CS,CSCL,CSGN,CSPN,CSR,CSS,CY,CZERO,C1,C2,DAI,PHI,PHID,RZ,
  !    *S1,S2,Y,Z,ZB,ZETA1,ZETA1D,ZETA2,ZETA2D,ZN,ZR
  REAL(8) :: aarg, aic, aii, air, Alim, ang, aphi, argdi, &
    argdr, argi, argr, asc, ascle, asumdi, asumdr, &
    asumi, asumr, bry, bsumdi, bsumdr, bsumi, bsumr, &
    car, cipi, cipr, cki, ckr, coner, crsc, cr1i, &
    cr1r, cr2i, cr2r, cscl, csgni, csi, cspni, cspnr, &
    csr, csrr, cssr, cyi, cyr, c1i, c1r, c2i, c2m, &
    c2r, daii, dair, Elim, fmr, fn, fnf, Fnu, hpi, &
    phidi, phidr, phii, phir, pi, pti, ptr, rast, &
    razr, rs1, rzi, rzr, sar, sgn, sti, str, s1i, &
    s1r, s2i, s2r, Tol, Yi, Yr, yy, zbi, zbr, &
    zeroi, zeror, zeta1i, zeta1r, zeta2i, zeta2r, &
    zet1di, zet1dr, zet2di, zet2dr, Zi, zni, znr, Zr, &
    zri, zrr, D1MACH, ZABS
  INTEGER i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, &
    kk, Kode, Mr, N, nai, ndai, nw, Nz, idum, j, ipard, ic
  DIMENSION bry(3), Yr(N), Yi(N), asumr(2), asumi(2), bsumr(2), &
    bsumi(2), phir(2), phii(2), argr(2), argi(2), zeta1r(2), &
    zeta1i(2), zeta2r(2), zeta2i(2), cyr(2), cyi(2), cipr(4), &
    cipi(4), cssr(3), csrr(3)
  EXTERNAL ZABS
  DATA zeror, zeroi, coner, cr1r, cr1i, cr2r, cr2i/0.0D0, 0.0D0, &
    1.0D0, 1.0D0, 1.73205080756887729D0, -0.5D0, &
    -8.66025403784438647D-01/
  DATA hpi, pi, aic/1.57079632679489662D+00, 3.14159265358979324D+00, &
    1.26551212348464539D+00/
  DATA cipr(1), cipi(1), cipr(2), cipi(2), cipr(3), cipi(3), cipr(4), &
    cipi(4)/1.0D0, 0.0D0, 0.0D0, -1.0D0, -1.0D0, 0.0D0, 0.0D0, &
    1.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZUNK2
  kdflg = 1
  Nz = 0
  !-----------------------------------------------------------------------
  !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
  !     THE UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  cscl = 1.0D0/Tol
  crsc = Tol
  cssr(1) = cscl
  cssr(2) = coner
  cssr(3) = crsc
  csrr(1) = crsc
  csrr(2) = coner
  csrr(3) = cscl
  bry(1) = 1.0D+3*D1MACH(1)/Tol
  bry(2) = 1.0D0/bry(1)
  bry(3) = D1MACH(2)
  zrr = Zr
  zri = Zi
  IF ( Zr<0.0D0 ) THEN
    zrr = -Zr
    zri = -Zi
  ENDIF
  yy = zri
  znr = zri
  zni = -zrr
  zbr = zrr
  zbi = zri
  inu = INT( Fnu )
  fnf = Fnu - inu
  ang = -hpi*fnf
  car = COS(ang)
  sar = SIN(ang)
  c2r = hpi*sar
  c2i = -hpi*car
  kk = MOD(inu,4) + 1
  str = c2r*cipr(kk) - c2i*cipi(kk)
  sti = c2r*cipi(kk) + c2i*cipr(kk)
  csr = cr1r*str - cr1i*sti
  csi = cr1r*sti + cr1i*str
  IF ( yy<=0.0D0 ) THEN
    znr = -znr
    zbi = -zbi
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
    CALL ZUNHJ(znr,zni,fn,0,Tol,phir(j),phii(j),argr(j),argi(j),zeta1r(j),&
      zeta1i(j),zeta2r(j),zeta2i(j),asumr(j),asumi(j),bsumr(j),&
      bsumi(j))
    IF ( Kode==1 ) THEN
      s1r = zeta1r(j) - zeta2r(j)
      s1i = zeta1i(j) - zeta2i(j)
    ELSE
      str = zbr + zeta2r(j)
      sti = zbi + zeta2i(j)
      rast = fn/ZABS(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = zeta1r(j) - str
      s1i = zeta1i(j) - sti
    ENDIF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = s1r
    IF ( ABS(rs1)<=Elim ) THEN
      IF ( kdflg==1 ) kflag = 2
      IF ( ABS(rs1)>=Alim ) THEN
        !-----------------------------------------------------------------------
        !     REFINE  TEST AND SCALE
        !-----------------------------------------------------------------------
        aphi = ZABS(phir(j),phii(j))
        aarg = ZABS(argr(j),argi(j))
        rs1 = rs1 + LOG(aphi) - 0.25D0*LOG(aarg) - aic
        IF ( ABS(rs1)>Elim ) GOTO 50
        IF ( kdflg==1 ) kflag = 1
        IF ( rs1>=0.0D0 ) THEN
          IF ( kdflg==1 ) kflag = 3
        ENDIF
      ENDIF
      !-----------------------------------------------------------------------
      !     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
      !     EXPONENT EXTREMES
      !-----------------------------------------------------------------------
      c2r = argr(j)*cr2r - argi(j)*cr2i
      c2i = argr(j)*cr2i + argi(j)*cr2r
      CALL ZAIRY(c2r,c2i,0,2,air,aii,nai,idum)
      CALL ZAIRY(c2r,c2i,1,2,dair,daii,ndai,idum)
      str = dair*bsumr(j) - daii*bsumi(j)
      sti = dair*bsumi(j) + daii*bsumr(j)
      ptr = str*cr2r - sti*cr2i
      pti = str*cr2i + sti*cr2r
      str = ptr + (air*asumr(j)-aii*asumi(j))
      sti = pti + (air*asumi(j)+aii*asumr(j))
      ptr = str*phir(j) - sti*phii(j)
      pti = str*phii(j) + sti*phir(j)
      s2r = ptr*csr - pti*csi
      s2i = ptr*csi + pti*csr
      str = EXP(s1r)*cssr(kflag)
      s1r = str*COS(s1i)
      s1i = str*SIN(s1i)
      str = s2r*s1r - s2i*s1i
      s2i = s1r*s2i + s2r*s1i
      s2r = str
      IF ( kflag==1 ) THEN
        CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
        IF ( nw/=0 ) GOTO 50
      ENDIF
      IF ( yy<=0.0D0 ) s2i = -s2i
      cyr(kdflg) = s2r
      cyi(kdflg) = s2i
      Yr(i) = s2r*csrr(kflag)
      Yi(i) = s2i*csrr(kflag)
      str = csi
      csi = -csr
      csr = str
      IF ( kdflg==2 ) GOTO 100
      kdflg = 2
      CYCLE
    ENDIF
    50     IF ( rs1>0.0D0 ) GOTO 600
    !-----------------------------------------------------------------------
    !     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
    !-----------------------------------------------------------------------
    IF ( Zr<0.0D0 ) GOTO 600
    kdflg = 1
    Yr(i) = zeror
    Yi(i) = zeroi
    Nz = Nz + 1
    str = csi
    csi = -csr
    csr = str
    IF ( i/=1 ) THEN
      IF ( (Yr(i-1)/=zeror).OR.(Yi(i-1)/=zeroi) ) THEN
        Yr(i-1) = zeror
        Yi(i-1) = zeroi
        Nz = Nz + 1
      ENDIF
    ENDIF
  ENDDO
  i = N
  100  razr = 1.0D0/ZABS(zrr,zri)
  str = zrr*razr
  sti = -zri*razr
  rzr = (str+str)*razr
  rzi = (sti+sti)*razr
  ckr = fn*rzr
  cki = fn*rzi
  ib = i + 1
  IF ( N<ib ) GOTO 300
  !-----------------------------------------------------------------------
  !     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
  !     ON UNDERFLOW.
  !-----------------------------------------------------------------------
  fn = Fnu + (N-1)
  ipard = 1
  IF ( Mr/=0 ) ipard = 0
  CALL ZUNHJ(znr,zni,fn,ipard,Tol,phidr,phidi,argdr,argdi,zet1dr,zet1di,&
    zet2dr,zet2di,asumdr,asumdi,bsumdr,bsumdi)
  IF ( Kode==1 ) THEN
    s1r = zet1dr - zet2dr
    s1i = zet1di - zet2di
  ELSE
    str = zbr + zet2dr
    sti = zbi + zet2di
    rast = fn/ZABS(str,sti)
    str = str*rast*rast
    sti = -sti*rast*rast
    s1r = zet1dr - str
    s1i = zet1di - sti
  ENDIF
  rs1 = s1r
  IF ( ABS(rs1)<=Elim ) THEN
    IF ( ABS(rs1)<Alim ) GOTO 200
    !-----------------------------------------------------------------------
    !     REFINE ESTIMATE AND TEST
    !-----------------------------------------------------------------------
    aphi = ZABS(phidr,phidi)
    rs1 = rs1 + LOG(aphi)
    IF ( ABS(rs1)<Elim ) GOTO 200
  ENDIF
  IF ( rs1>0.0D0 ) GOTO 600
  !-----------------------------------------------------------------------
  !     FOR ZR.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
  !-----------------------------------------------------------------------
  IF ( Zr<0.0D0 ) GOTO 600
  Nz = N
  DO i = 1, N
    Yr(i) = zeror
    Yi(i) = zeroi
  ENDDO
  RETURN
  200  s1r = cyr(1)
  s1i = cyi(1)
  s2r = cyr(2)
  s2i = cyi(2)
  c1r = csrr(kflag)
  ascle = bry(kflag)
  DO i = ib, N
    c2r = s2r
    c2i = s2i
    s2r = ckr*c2r - cki*c2i + s1r
    s2i = ckr*c2i + cki*c2r + s1i
    s1r = c2r
    s1i = c2i
    ckr = ckr + rzr
    cki = cki + rzi
    c2r = s2r*c1r
    c2i = s2i*c1r
    Yr(i) = c2r
    Yi(i) = c2i
    IF ( kflag<3 ) THEN
      str = ABS(c2r)
      sti = ABS(c2i)
      c2m = MAX(str,sti)
      IF ( c2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*c1r
        s1i = s1i*c1r
        s2r = c2r
        s2i = c2i
        s1r = s1r*cssr(kflag)
        s1i = s1i*cssr(kflag)
        s2r = s2r*cssr(kflag)
        s2i = s2i*cssr(kflag)
        c1r = csrr(kflag)
      ENDIF
    ENDIF
  ENDDO
  300 CONTINUE
  IF ( Mr==0 ) RETURN
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION FOR RE(Z).LT.0.0D0
  !-----------------------------------------------------------------------
  Nz = 0
  fmr = Mr
  sgn = -DSIGN(pi,fmr)
  !-----------------------------------------------------------------------
  !     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
  !-----------------------------------------------------------------------
  csgni = sgn
  IF ( yy<=0.0D0 ) csgni = -csgni
  ifn = inu + N - 1
  ang = fnf*sgn
  cspnr = COS(ang)
  cspni = SIN(ang)
  IF ( MOD(ifn,2)/=0 ) THEN
    cspnr = -cspnr
    cspni = -cspni
  ENDIF
  !-----------------------------------------------------------------------
  !     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION. I(FNU,Z) IS
  !     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
  !     QUADRANT. FOURTH QUADRANT VALUES (YY.LE.0.0E0) ARE COMPUTED BY
  !     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
  !-----------------------------------------------------------------------
  csr = sar*csgni
  csi = car*csgni
  in = MOD(ifn,4) + 1
  c2r = cipr(in)
  c2i = cipi(in)
  str = csr*c2r + csi*c2i
  csi = -csr*c2i + csi*c2r
  csr = str
  asc = bry(1)
  iuf = 0
  kk = N
  kdflg = 1
  ib = ib - 1
  ic = ib - 1
  DO k = 1, N
    fn = Fnu + (kk-1)
    !-----------------------------------------------------------------------
    !     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
    !     FUNCTION ABOVE
    !-----------------------------------------------------------------------
    IF ( N>2 ) THEN
      IF ( (kk==N).AND.(ib<N) ) GOTO 350
      IF ( (kk/=ib).AND.(kk/=ic) ) THEN
        CALL ZUNHJ(znr,zni,fn,0,Tol,phidr,phidi,argdr,argdi,zet1dr,zet1di,&
          zet2dr,zet2di,asumdr,asumdi,bsumdr,bsumdi)
        GOTO 350
      ENDIF
    ENDIF
    phidr = phir(j)
    phidi = phii(j)
    argdr = argr(j)
    argdi = argi(j)
    zet1dr = zeta1r(j)
    zet1di = zeta1i(j)
    zet2dr = zeta2r(j)
    zet2di = zeta2i(j)
    asumdr = asumr(j)
    asumdi = asumi(j)
    bsumdr = bsumr(j)
    bsumdi = bsumi(j)
    j = 3 - j
    350    IF ( Kode==1 ) THEN
    s1r = -zet1dr + zet2dr
    s1i = -zet1di + zet2di
  ELSE
    str = zbr + zet2dr
    sti = zbi + zet2di
    rast = fn/ZABS(str,sti)
    str = str*rast*rast
    sti = -sti*rast*rast
    s1r = -zet1dr + str
    s1i = -zet1di + sti
  ENDIF
  !-----------------------------------------------------------------------
  !     TEST FOR UNDERFLOW AND OVERFLOW
  !-----------------------------------------------------------------------
  rs1 = s1r
  IF ( ABS(rs1)>Elim ) GOTO 450
  IF ( kdflg==1 ) iflag = 2
  IF ( ABS(rs1)>=Alim ) THEN
    !-----------------------------------------------------------------------
    !     REFINE  TEST AND SCALE
    !-----------------------------------------------------------------------
    aphi = ZABS(phidr,phidi)
    aarg = ZABS(argdr,argdi)
    rs1 = rs1 + LOG(aphi) - 0.25D0*LOG(aarg) - aic
    IF ( ABS(rs1)>Elim ) GOTO 450
    IF ( kdflg==1 ) iflag = 1
    IF ( rs1>=0.0D0 ) THEN
      IF ( kdflg==1 ) iflag = 3
    ENDIF
  ENDIF
  CALL ZAIRY(argdr,argdi,0,2,air,aii,nai,idum)
  CALL ZAIRY(argdr,argdi,1,2,dair,daii,ndai,idum)
  str = dair*bsumdr - daii*bsumdi
  sti = dair*bsumdi + daii*bsumdr
  str = str + (air*asumdr-aii*asumdi)
  sti = sti + (air*asumdi+aii*asumdr)
  ptr = str*phidr - sti*phidi
  pti = str*phidi + sti*phidr
  s2r = ptr*csr - pti*csi
  s2i = ptr*csi + pti*csr
  str = EXP(s1r)*cssr(iflag)
  s1r = str*COS(s1i)
  s1i = str*SIN(s1i)
  str = s2r*s1r - s2i*s1i
  s2i = s2r*s1i + s2i*s1r
  s2r = str
  IF ( iflag==1 ) THEN
    CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
    IF ( nw/=0 ) THEN
      s2r = zeror
      s2i = zeroi
    ENDIF
  ENDIF
  400    IF ( yy<=0.0D0 ) s2i = -s2i
  cyr(kdflg) = s2r
  cyi(kdflg) = s2i
  c2r = s2r
  c2i = s2i
  s2r = s2r*csrr(iflag)
  s2i = s2i*csrr(iflag)
  !-----------------------------------------------------------------------
  !     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
  !-----------------------------------------------------------------------
  s1r = Yr(kk)
  s1i = Yi(kk)
  IF ( Kode/=1 ) THEN
    CALL ZS1S2(zrr,zri,s1r,s1i,s2r,s2i,nw,asc,Alim,iuf)
    Nz = Nz + nw
  ENDIF
  Yr(kk) = s1r*cspnr - s1i*cspni + s2r
  Yi(kk) = s1r*cspni + s1i*cspnr + s2i
  kk = kk - 1
  cspnr = -cspnr
  cspni = -cspni
  str = csi
  csi = -csr
  csr = str
  IF ( c2r/=0.0D0.OR.c2i/=0.0D0 ) THEN
    IF ( kdflg==2 ) GOTO 500
    kdflg = 2
    CYCLE
  ELSE
    kdflg = 1
    CYCLE
  ENDIF
  450    IF ( rs1>0.0D0 ) GOTO 600
  s2r = zeror
  s2i = zeroi
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
  s1r = cyr(1)
  s1i = cyi(1)
  s2r = cyr(2)
  s2i = cyi(2)
  csr = csrr(iflag)
  ascle = bry(iflag)
  fn = inu + il
  DO i = 1, il
    c2r = s2r
    c2i = s2i
    s2r = s1r + (fn+fnf)*(rzr*c2r-rzi*c2i)
    s2i = s1i + (fn+fnf)*(rzr*c2i+rzi*c2r)
    s1r = c2r
    s1i = c2i
    fn = fn - 1.0D0
    c2r = s2r*csr
    c2i = s2i*csr
    ckr = c2r
    cki = c2i
    c1r = Yr(kk)
    c1i = Yi(kk)
    IF ( Kode/=1 ) THEN
      CALL ZS1S2(zrr,zri,c1r,c1i,c2r,c2i,nw,asc,Alim,iuf)
      Nz = Nz + nw
    ENDIF
    Yr(kk) = c1r*cspnr - c1i*cspni + c2r
    Yi(kk) = c1r*cspni + c1i*cspnr + c2i
    kk = kk - 1
    cspnr = -cspnr
    cspni = -cspni
    IF ( iflag<3 ) THEN
      c2r = ABS(ckr)
      c2i = ABS(cki)
      c2m = MAX(c2r,c2i)
      IF ( c2m>ascle ) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1r = s1r*csr
        s1i = s1i*csr
        s2r = ckr
        s2i = cki
        s1r = s1r*cssr(iflag)
        s1i = s1i*cssr(iflag)
        s2r = s2r*cssr(iflag)
        s2i = s2i*cssr(iflag)
        csr = csrr(iflag)
      ENDIF
    ENDIF
  ENDDO
  RETURN
  600  Nz = -1
END SUBROUTINE ZUNK2