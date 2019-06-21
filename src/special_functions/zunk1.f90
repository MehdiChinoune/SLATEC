!** ZUNK1
SUBROUTINE ZUNK1(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUNK1-A, ZUNK1-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
  !     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
  !     UNIFORM ASYMPTOTIC EXPANSION.
  !     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
  !     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
  !
  !***
  ! **See also:**  ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZS1S2, ZUCHK, ZUNIK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  !     COMPLEX CFN,CK,CONE,CRSC,CS,CSCL,CSGN,CSPN,CSR,CSS,CWRK,CY,CZERO,
  !    *C1,C2,PHI,PHID,RZ,SUM,SUMD,S1,S2,Y,Z,ZETA1,ZETA1D,ZETA2,ZETA2D,ZR
  INTEGER :: i, ib, iflag, ifn, il, init(2), inu, iuf, k, kdflg, kflag, &
    kk, Kode, Mr, N, nw, Nz, initd, ic, ipard, j, m
  REAL(DP) :: Alim, ang, aphi, asc, ascle, bry(3), cki, ckr, crsc, cscl, csgni, &
    cspni, cspnr, csr, csrr(3), cssr(3), cwrki(16,3), cwrkr(16,3), cyi(2), cyr(2), &
    c1i, c1r, c2i, c2m, c2r, Elim, fmr, fn, fnf, Fnu, phidi, phidr, phii(2), &
    phir(2), rast, razr, rs1, rzi, rzr, sgn, sti, str, sumdi, sumdr, sumi(2), &
    sumr(2), s1i, s1r, s2i, s2r, Tol, Yi(N), Yr(N), zeta1i(2), zeta1r(2), &
    zeta2i(2), zeta2r(2), zet1di, zet1dr, zet2di, zet2dr, Zi, Zr, zri, zrr
  REAL(DP), PARAMETER :: zeror = 0._DP, zeroi = 0._DP, coner = 1._DP
  REAL(DP), PARAMETER :: pi = 3.14159265358979324_DP
  !* FIRST EXECUTABLE STATEMENT  ZUNK1
  kdflg = 1
  Nz = 0
  !-----------------------------------------------------------------------
  !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
  !     THE UNDERFLOW LIMIT
  !-----------------------------------------------------------------------
  cscl = 1._DP/Tol
  crsc = Tol
  cssr(1) = cscl
  cssr(2) = coner
  cssr(3) = crsc
  csrr(1) = crsc
  csrr(2) = coner
  csrr(3) = cscl
  bry(1) = 1.E3_DP*D1MACH(1)/Tol
  bry(2) = 1._DP/bry(1)
  bry(3) = D1MACH(2)
  zrr = Zr
  zri = Zi
  IF( Zr<0._DP ) THEN
    zrr = -Zr
    zri = -Zi
  END IF
  j = 2
  DO i = 1, N
    !-----------------------------------------------------------------------
    !     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
    !-----------------------------------------------------------------------
    j = 3 - j
    fn = Fnu + (i-1)
    init(j) = 0
    CALL ZUNIK(zrr,zri,fn,2,0,Tol,init(j),phir(j),phii(j),zeta1r(j),&
      zeta1i(j),zeta2r(j),zeta2i(j),sumr(j),sumi(j),cwrkr(1,j),cwrki(1,j))
    IF( Kode==1 ) THEN
      s1r = zeta1r(j) - zeta2r(j)
      s1i = zeta1i(j) - zeta2i(j)
    ELSE
      str = zrr + zeta2r(j)
      sti = zri + zeta2i(j)
      rast = fn/ZABS(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = zeta1r(j) - str
      s1i = zeta1i(j) - sti
    END IF
    rs1 = s1r
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    IF( ABS(rs1)<=Elim ) THEN
      IF( kdflg==1 ) kflag = 2
      IF( ABS(rs1)>=Alim ) THEN
        !-----------------------------------------------------------------------
        !     REFINE  TEST AND SCALE
        !-----------------------------------------------------------------------
        aphi = ZABS(phir(j),phii(j))
        rs1 = rs1 + LOG(aphi)
        IF( ABS(rs1)>Elim ) GOTO 50
        IF( kdflg==1 ) kflag = 1
        IF( rs1>=0._DP ) THEN
          IF( kdflg==1 ) kflag = 3
        END IF
      END IF
      !-----------------------------------------------------------------------
      !     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
      !     EXPONENT EXTREMES
      !-----------------------------------------------------------------------
      s2r = phir(j)*sumr(j) - phii(j)*sumi(j)
      s2i = phir(j)*sumi(j) + phii(j)*sumr(j)
      str = EXP(s1r)*cssr(kflag)
      s1r = str*COS(s1i)
      s1i = str*SIN(s1i)
      str = s2r*s1r - s2i*s1i
      s2i = s1r*s2i + s2r*s1i
      s2r = str
      IF( kflag==1 ) THEN
        CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
        IF( nw/=0 ) GOTO 50
      END IF
      cyr(kdflg) = s2r
      cyi(kdflg) = s2i
      Yr(i) = s2r*csrr(kflag)
      Yi(i) = s2i*csrr(kflag)
      IF( kdflg==2 ) GOTO 100
      kdflg = 2
      CYCLE
    END IF
    50  IF( rs1>0._DP ) GOTO 600
    !-----------------------------------------------------------------------
    !     FOR ZR<0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
    !-----------------------------------------------------------------------
    IF( Zr<0._DP ) GOTO 600
    kdflg = 1
    Yr(i) = zeror
    Yi(i) = zeroi
    Nz = Nz + 1
    IF( i/=1 ) THEN
      IF( (Yr(i-1)/=zeror) .OR. (Yi(i-1)/=zeroi) ) THEN
        Yr(i-1) = zeror
        Yi(i-1) = zeroi
        Nz = Nz + 1
      END IF
    END IF
  END DO
  i = N
  100  razr = 1._DP/ZABS(zrr,zri)
  str = zrr*razr
  sti = -zri*razr
  rzr = (str+str)*razr
  rzi = (sti+sti)*razr
  ckr = fn*rzr
  cki = fn*rzi
  ib = i + 1
  IF( N<ib ) GOTO 300
  !-----------------------------------------------------------------------
  !     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW. SET SEQUENCE TO ZERO
  !     ON UNDERFLOW.
  !-----------------------------------------------------------------------
  fn = Fnu + (N-1)
  ipard = 1
  IF( Mr/=0 ) ipard = 0
  initd = 0
  CALL ZUNIK(zrr,zri,fn,2,ipard,Tol,initd,phidr,phidi,zet1dr,zet1di,zet2dr,&
    zet2di,sumdr,sumdi,cwrkr(1,3),cwrki(1,3))
  IF( Kode==1 ) THEN
    s1r = zet1dr - zet2dr
    s1i = zet1di - zet2di
  ELSE
    str = zrr + zet2dr
    sti = zri + zet2di
    rast = fn/ZABS(str,sti)
    str = str*rast*rast
    sti = -sti*rast*rast
    s1r = zet1dr - str
    s1i = zet1di - sti
  END IF
  rs1 = s1r
  IF( ABS(rs1)<=Elim ) THEN
    IF( ABS(rs1)<Alim ) GOTO 200
    !-----------------------------------------------------------------------
    !     REFINE ESTIMATE AND TEST
    !-----------------------------------------------------------------------
    aphi = ZABS(phidr,phidi)
    rs1 = rs1 + LOG(aphi)
    IF( ABS(rs1)<Elim ) GOTO 200
  END IF
  IF( ABS(rs1)>0._DP ) GOTO 600
  !-----------------------------------------------------------------------
  !     FOR ZR<0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
  !-----------------------------------------------------------------------
  IF( Zr<0._DP ) GOTO 600
  Nz = N
  DO i = 1, N
    Yr(i) = zeror
    Yi(i) = zeroi
  END DO
  RETURN
  !-----------------------------------------------------------------------
  !     FORWARD RECUR FOR REMAINDER OF THE SEQUENCE
  !-----------------------------------------------------------------------
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
    IF( kflag<3 ) THEN
      str = ABS(c2r)
      sti = ABS(c2i)
      c2m = MAX(str,sti)
      IF( c2m>ascle ) THEN
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
      END IF
    END IF
  END DO
  300 CONTINUE
  IF( Mr==0 ) RETURN
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION FOR RE(Z)<0._DP
  !-----------------------------------------------------------------------
  Nz = 0
  fmr = Mr
  sgn = -SIGN(pi,fmr)
  !-----------------------------------------------------------------------
  !     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
  !-----------------------------------------------------------------------
  csgni = sgn
  inu = INT( Fnu )
  fnf = Fnu - inu
  ifn = inu + N - 1
  ang = fnf*sgn
  cspnr = COS(ang)
  cspni = SIN(ang)
  IF( MOD(ifn,2)/=0 ) THEN
    cspnr = -cspnr
    cspni = -cspni
  END IF
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
    m = 3
    IF( N>2 ) THEN
      IF( (kk==N) .AND. (ib<N) ) GOTO 350
      IF( (kk/=ib) .AND. (kk/=ic) ) THEN
        initd = 0
        GOTO 350
      END IF
    END IF
    initd = init(j)
    phidr = phir(j)
    phidi = phii(j)
    zet1dr = zeta1r(j)
    zet1di = zeta1i(j)
    zet2dr = zeta2r(j)
    zet2di = zeta2i(j)
    sumdr = sumr(j)
    sumdi = sumi(j)
    m = j
    j = 3 - j
    350  CALL ZUNIK(zrr,zri,fn,1,0,Tol,initd,phidr,phidi,zet1dr,zet1di,zet2dr,&
      zet2di,sumdr,sumdi,cwrkr(1,m),cwrki(1,m))
    IF( Kode==1 ) THEN
      s1r = -zet1dr + zet2dr
      s1i = -zet1di + zet2di
    ELSE
      str = zrr + zet2dr
      sti = zri + zet2di
      rast = fn/ZABS(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = -zet1dr + str
      s1i = -zet1di + sti
    END IF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = s1r
    IF( ABS(rs1)>Elim ) GOTO 450
    IF( kdflg==1 ) iflag = 2
    IF( ABS(rs1)>=Alim ) THEN
      !-----------------------------------------------------------------------
      !     REFINE  TEST AND SCALE
      !-----------------------------------------------------------------------
      aphi = ZABS(phidr,phidi)
      rs1 = rs1 + LOG(aphi)
      IF( ABS(rs1)>Elim ) GOTO 450
      IF( kdflg==1 ) iflag = 1
      IF( rs1>=0._DP ) THEN
        IF( kdflg==1 ) iflag = 3
      END IF
    END IF
    str = phidr*sumdr - phidi*sumdi
    sti = phidr*sumdi + phidi*sumdr
    s2r = -csgni*sti
    s2i = csgni*str
    str = EXP(s1r)*cssr(iflag)
    s1r = str*COS(s1i)
    s1i = str*SIN(s1i)
    str = s2r*s1r - s2i*s1i
    s2i = s2r*s1i + s2i*s1r
    s2r = str
    IF( iflag==1 ) THEN
      CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
      IF( nw/=0 ) THEN
        s2r = zeror
        s2i = zeroi
      END IF
    END IF
    400  cyr(kdflg) = s2r
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
    IF( Kode/=1 ) THEN
      CALL ZS1S2(zrr,zri,s1r,s1i,s2r,s2i,nw,asc,Alim,iuf)
      Nz = Nz + nw
    END IF
    Yr(kk) = s1r*cspnr - s1i*cspni + s2r
    Yi(kk) = cspnr*s1i + cspni*s1r + s2i
    kk = kk - 1
    cspnr = -cspnr
    cspni = -cspni
    IF( c2r/=0._DP .OR. c2i/=0._DP ) THEN
      IF( kdflg==2 ) GOTO 500
      kdflg = 2
      CYCLE
    ELSE
      kdflg = 1
      CYCLE
    END IF
    450  IF( rs1>0._DP ) GOTO 600
    s2r = zeror
    s2i = zeroi
    GOTO 400
  END DO
  k = N
  500  il = N - k
  IF( il==0 ) RETURN
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
    fn = fn - 1._DP
    c2r = s2r*csr
    c2i = s2i*csr
    ckr = c2r
    cki = c2i
    c1r = Yr(kk)
    c1i = Yi(kk)
    IF( Kode/=1 ) THEN
      CALL ZS1S2(zrr,zri,c1r,c1i,c2r,c2i,nw,asc,Alim,iuf)
      Nz = Nz + nw
    END IF
    Yr(kk) = c1r*cspnr - c1i*cspni + c2r
    Yi(kk) = c1r*cspni + c1i*cspnr + c2i
    kk = kk - 1
    cspnr = -cspnr
    cspni = -cspni
    IF( iflag<3 ) THEN
      c2r = ABS(ckr)
      c2i = ABS(cki)
      c2m = MAX(c2r,c2i)
      IF( c2m>ascle ) THEN
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
      END IF
    END IF
  END DO
  RETURN
  600  Nz = -1
END SUBROUTINE ZUNK1
