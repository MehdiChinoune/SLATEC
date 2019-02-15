!DECK CBKNU
SUBROUTINE CBKNU(Z,Fnu,Kode,N,Y,Nz,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CBKNU
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CAIRY, CBESH, CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CBKNU-A, ZBKNU-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
  !
  !***SEE ALSO  CAIRY, CBESH, CBESI, CBESK
  !***ROUTINES CALLED  CKSCL, CSHCH, CUCHK, GAMLN, I1MACH, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CBKNU
  !
  COMPLEX cch, ck, coef, cone, crsc, cs, cscl, csh, csr, css, &
    ctwo, cz, czero, f, fmu, p, pt, p1, p2, q, rz, smu, &
    st, s1, s2, Y, Z, zd, celm, cy
  REAL aa, ak, Alim, ascle, a1, a2, bb, bk, bry, caz, cc, dnu, &
    dnu2, Elim, etest, fc, fhs, fk, fks, Fnu, fpi, g1, g2, &
    hpi, pi, p2i, p2m, p2r, rk, rthpi, r1, s, spi, tm, Tol, &
    tth, t1, t2, xx, yy, GAMLN, R1MACH, helim, elm, xd, yd, &
    alas, as
  INTEGER i, idum, iflag, inu, k, kflag, kk, kmax, Kode, koded, &
    N, Nz, I1MACH, nw, j, ic, inub
  DIMENSION bry(3), cc(8), css(3), csr(3), Y(N), cy(2)
  !
  DATA kmax/30/
  DATA r1/2.0E0/
  DATA czero, cone, ctwo/(0.0E0,0.0E0), (1.0E0,0.0E0), (2.0E0,0.0E0)/
  !
  DATA pi, rthpi, spi, hpi, fpi, tth/3.14159265358979324E0, &
    1.25331413731550025E0, 1.90985931710274403E0, &
    1.57079632679489662E0, 1.89769999331517738E0, &
    6.66666666666666666E-01/
  !
  DATA cc(1), cc(2), cc(3), cc(4), cc(5), cc(6), cc(7), &
    cc(8)/5.77215664901532861E-01, -4.20026350340952355E-02, &
    -4.21977345555443367E-02, 7.21894324666309954E-03, &
    -2.15241674114950973E-04, -2.01348547807882387E-05, &
    1.13302723198169588E-06, 6.11609510448141582E-09/
  !
  !***FIRST EXECUTABLE STATEMENT  CBKNU
  xx = REAL(Z)
  yy = AIMAG(Z)
  caz = ABS(Z)
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
  Nz = 0
  iflag = 0
  koded = Kode
  rz = ctwo/Z
  inu = Fnu + 0.5E0
  dnu = Fnu - inu
  IF ( ABS(dnu)/=0.5E0 ) THEN
    dnu2 = 0.0E0
    IF ( ABS(dnu)>Tol ) dnu2 = dnu*dnu
    IF ( caz<=r1 ) THEN
      !-----------------------------------------------------------------------
      !     SERIES FOR ABS(Z).LE.R1
      !-----------------------------------------------------------------------
      fc = 1.0E0
      smu = CLOG(rz)
      fmu = smu*CMPLX(dnu,0.0E0)
      CALL CSHCH(fmu,csh,cch)
      IF ( dnu/=0.0E0 ) THEN
        fc = dnu*pi
        fc = fc/SIN(fc)
        smu = csh*CMPLX(1.0E0/dnu,0.0E0)
      ENDIF
      a2 = 1.0E0 + dnu
      !-----------------------------------------------------------------------
      !     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
      !-----------------------------------------------------------------------
      t2 = EXP(-GAMLN(a2,idum))
      t1 = 1.0E0/(t2*fc)
      IF ( ABS(dnu)>0.1E0 ) THEN
        g1 = (t1-t2)/(dnu+dnu)
      ELSE
        !-----------------------------------------------------------------------
        !     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
        !-----------------------------------------------------------------------
        ak = 1.0E0
        s = cc(1)
        DO k = 2, 8
          ak = ak*dnu2
          tm = cc(k)*ak
          s = s + tm
          IF ( ABS(tm)<Tol ) EXIT
        ENDDO
        g1 = -s
      ENDIF
      g2 = 0.5E0*(t1+t2)*fc
      g1 = g1*fc
      f = CMPLX(g1,0.0E0)*cch + smu*CMPLX(g2,0.0E0)
      pt = CEXP(fmu)
      p = CMPLX(0.5E0/t2,0.0E0)*pt
      q = CMPLX(0.5E0/t1,0.0E0)/pt
      s1 = f
      s2 = p
      ak = 1.0E0
      a1 = 1.0E0
      ck = cone
      bk = 1.0E0 - dnu2
      IF ( inu>0.OR.N>1 ) THEN
        !-----------------------------------------------------------------------
        !     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
        !-----------------------------------------------------------------------
        IF ( caz>=Tol ) THEN
          cz = Z*Z*CMPLX(0.25E0,0.0E0)
          t1 = 0.25E0*caz*caz
          DO
            f = (f*CMPLX(ak,0.0E0)+p+q)*CMPLX(1.0E0/bk,0.0E0)
            p = p*CMPLX(1.0E0/(ak-dnu),0.0E0)
            q = q*CMPLX(1.0E0/(ak+dnu),0.0E0)
            rk = 1.0E0/ak
            ck = ck*cz*CMPLX(rk,0.0E0)
            s1 = s1 + ck*f
            s2 = s2 + ck*(p-f*CMPLX(ak,0.0E0))
            a1 = a1*t1*rk
            bk = bk + ak + ak + 1.0E0
            ak = ak + 1.0E0
            IF ( a1<=Tol ) EXIT
          ENDDO
        ENDIF
        kflag = 2
        bk = REAL(smu)
        a1 = Fnu + 1.0E0
        ak = a1*ABS(bk)
        IF ( ak>Alim ) kflag = 3
        p2 = s2*css(kflag)
        s2 = p2*rz
        s1 = s1*css(kflag)
        IF ( koded/=1 ) THEN
          f = CEXP(Z)
          s1 = s1*f
          s2 = s2*f
        ENDIF
        GOTO 200
      ELSE
        !-----------------------------------------------------------------------
        !     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
        !-----------------------------------------------------------------------
        IF ( caz>=Tol ) THEN
          cz = Z*Z*CMPLX(0.25E0,0.0E0)
          t1 = 0.25E0*caz*caz
          DO
            f = (f*CMPLX(ak,0.0E0)+p+q)*CMPLX(1.0E0/bk,0.0E0)
            p = p*CMPLX(1.0E0/(ak-dnu),0.0E0)
            q = q*CMPLX(1.0E0/(ak+dnu),0.0E0)
            rk = 1.0E0/ak
            ck = ck*cz*CMPLX(rk,0.0)
            s1 = s1 + ck*f
            a1 = a1*t1*rk
            bk = bk + ak + ak + 1.0E0
            ak = ak + 1.0E0
            IF ( a1<=Tol ) EXIT
          ENDDO
        ENDIF
        Y(1) = s1
        IF ( koded==1 ) RETURN
        Y(1) = s1*CEXP(Z)
        RETURN
      ENDIF
    ENDIF
  ENDIF
  !-----------------------------------------------------------------------
  !     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
  !     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
  !     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
  !     RECURSION
  !-----------------------------------------------------------------------
  coef = CMPLX(rthpi,0.0E0)/CSQRT(Z)
  kflag = 2
  IF ( koded/=2 ) THEN
    IF ( xx>Alim ) THEN
      !-----------------------------------------------------------------------
      !     SCALE BY EXP(Z), IFLAG = 1 CASES
      !-----------------------------------------------------------------------
      koded = 2
      iflag = 1
      kflag = 2
    ELSE
      !     BLANK LINE
      a1 = EXP(-xx)*REAL(css(kflag))
      pt = CMPLX(a1,0.0E0)*CMPLX(COS(yy),-SIN(yy))
      coef = coef*pt
    ENDIF
  ENDIF
  IF ( ABS(dnu)==0.5E0 ) GOTO 800
  !-----------------------------------------------------------------------
  !     MILLER ALGORITHM FOR ABS(Z).GT.R1
  !-----------------------------------------------------------------------
  ak = COS(pi*dnu)
  ak = ABS(ak)
  IF ( ak==0.0E0 ) GOTO 800
  fhs = ABS(0.25E0-dnu2)
  IF ( fhs==0.0E0 ) GOTO 800
  !-----------------------------------------------------------------------
  !     COMPUTE R2=F(E). IF ABS(Z).GE.R2, USE FORWARD RECURRENCE TO
  !     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
  !     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(11))=
  !     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
  !-----------------------------------------------------------------------
  t1 = (I1MACH(11)-1)*R1MACH(5)*3.321928094E0
  t1 = MAX(t1,12.0E0)
  t1 = MIN(t1,60.0E0)
  t2 = tth*t1 - 6.0E0
  IF ( xx/=0.0E0 ) THEN
    t1 = ATAN(yy/xx)
    t1 = ABS(t1)
  ELSE
    t1 = hpi
  ENDIF
  IF ( t2>caz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE BACKWARD INDEX K FOR ABS(Z).LT.R2
    !-----------------------------------------------------------------------
    a2 = SQRT(caz)
    ak = fpi*ak/(Tol*SQRT(a2))
    aa = 3.0E0*t1/(1.0E0+caz)
    bb = 14.7E0*t1/(28.0E0+caz)
    ak = (ALOG(ak)+caz*COS(aa)/(1.0E0+0.008E0*caz))/COS(bb)
    fk = 0.12125E0*ak*ak/caz + 1.5E0
  ELSE
    !-----------------------------------------------------------------------
    !     FORWARD RECURRENCE LOOP WHEN ABS(Z).GE.R2
    !-----------------------------------------------------------------------
    etest = ak/(pi*caz*Tol)
    fk = 1.0E0
    IF ( etest<1.0E0 ) GOTO 100
    fks = 2.0E0
    rk = caz + caz + 2.0E0
    a1 = 0.0E0
    a2 = 1.0E0
    DO i = 1, kmax
      ak = fhs/fks
      bk = rk/(fk+1.0E0)
      tm = a2
      a2 = bk*a2 - ak*a1
      a1 = tm
      rk = rk + 2.0E0
      fks = fks + fk + fk + 2.0E0
      fhs = fhs + fk + fk
      fk = fk + 1.0E0
      tm = ABS(a2)*fk
      IF ( etest<tm ) GOTO 50
    ENDDO
    Nz = -2
    GOTO 99999
    50     fk = fk + spi*t1*SQRT(t2/caz)
    fhs = ABS(0.25E0-dnu2)
  ENDIF
  100  k = fk
  !-----------------------------------------------------------------------
  !     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  fk = k
  fks = fk*fk
  p1 = czero
  p2 = CMPLX(Tol,0.0E0)
  cs = p2
  DO i = 1, k
    a1 = fks - fk
    a2 = (fks+fk)/(a1+fhs)
    rk = 2.0E0/(fk+1.0E0)
    t1 = (fk+xx)*rk
    t2 = yy*rk
    pt = p2
    p2 = (p2*CMPLX(t1,t2)-p1)*CMPLX(a2,0.0E0)
    p1 = pt
    cs = cs + p2
    fks = a1 - fk + 1.0E0
    fk = fk - 1.0E0
  ENDDO
  !-----------------------------------------------------------------------
  !     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
  !     SCALING
  !-----------------------------------------------------------------------
  tm = ABS(cs)
  pt = CMPLX(1.0E0/tm,0.0E0)
  s1 = pt*p2
  cs = CONJG(cs)*pt
  s1 = coef*s1*cs
  IF ( inu>0.OR.N>1 ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
    !-----------------------------------------------------------------------
    tm = ABS(p2)
    pt = CMPLX(1.0E0/tm,0.0E0)
    p1 = pt*p1
    p2 = CONJG(p2)*pt
    pt = p1*p2
    s2 = s1*(cone+(CMPLX(dnu+0.5E0,0.0E0)-pt)/Z)
  ELSE
    zd = Z
    IF ( iflag/=1 ) GOTO 400
    GOTO 700
  ENDIF
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
  !     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
  !-----------------------------------------------------------------------
  200  ck = CMPLX(dnu+1.0E0,0.0E0)*rz
  IF ( N==1 ) inu = inu - 1
  IF ( inu>0 ) THEN
    inub = 1
    IF ( iflag==1 ) THEN
      !-----------------------------------------------------------------------
      !     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
      !-----------------------------------------------------------------------
      helim = 0.5E0*Elim
      elm = EXP(-Elim)
      celm = CMPLX(elm,0.0)
      ascle = bry(1)
      zd = Z
      xd = xx
      yd = yy
      ic = -1
      j = 2
      DO i = 1, inu
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rz
        as = ABS(s2)
        alas = ALOG(as)
        p2r = -xd + alas
        IF ( p2r>=(-Elim) ) THEN
          p2 = -zd + CLOG(s2)
          p2r = REAL(p2)
          p2i = AIMAG(p2)
          p2m = EXP(p2r)/Tol
          p1 = CMPLX(p2m,0.0E0)*CMPLX(COS(p2i),SIN(p2i))
          CALL CUCHK(p1,nw,ascle,Tol)
          IF ( nw==0 ) THEN
            j = 3 - j
            cy(j) = p1
            IF ( ic==(i-1) ) GOTO 600
            ic = i
            CYCLE
          ENDIF
        ENDIF
        IF ( alas>=helim ) THEN
          xd = xd - Elim
          s1 = s1*celm
          s2 = s2*celm
          zd = CMPLX(xd,yd)
        ENDIF
      ENDDO
      IF ( N==1 ) s1 = s2
      GOTO 700
    ENDIF
  ELSE
    IF ( N==1 ) s1 = s2
    zd = Z
    IF ( iflag/=1 ) GOTO 400
    GOTO 700
  ENDIF
  300  p1 = csr(kflag)
  ascle = bry(kflag)
  DO i = inub, inu
    st = s2
    s2 = ck*s2 + s1
    s1 = st
    ck = ck + rz
    IF ( kflag<3 ) THEN
      p2 = s2*p1
      p2r = REAL(p2)
      p2i = AIMAG(p2)
      p2r = ABS(p2r)
      p2i = ABS(p2i)
      p2m = MAX(p2r,p2i)
      IF ( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
      ENDIF
    ENDIF
  ENDDO
  IF ( N==1 ) s1 = s2
  400  Y(1) = s1*csr(kflag)
  IF ( N==1 ) RETURN
  Y(2) = s2*csr(kflag)
  IF ( N==2 ) RETURN
  kk = 2
  500  kk = kk + 1
  IF ( kk>N ) RETURN
  p1 = csr(kflag)
  ascle = bry(kflag)
  DO i = kk, N
    p2 = s2
    s2 = ck*s2 + s1
    s1 = p2
    ck = ck + rz
    p2 = s2*p1
    Y(i) = p2
    IF ( kflag<3 ) THEN
      p2r = REAL(p2)
      p2i = AIMAG(p2)
      p2r = ABS(p2r)
      p2i = ABS(p2i)
      p2m = MAX(p2r,p2i)
      IF ( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
      ENDIF
    ENDIF
  ENDDO
  RETURN
  600  kflag = 1
  inub = i + 1
  s2 = cy(j)
  j = 3 - j
  s1 = cy(j)
  IF ( inub<=inu ) GOTO 300
  IF ( N==1 ) s1 = s2
  GOTO 400
  700  Y(1) = s1
  IF ( N/=1 ) Y(2) = s2
  ascle = bry(1)
  CALL CKSCL(zd,Fnu,N,Y,Nz,rz,ascle,Tol,Elim)
  inu = N - Nz
  IF ( inu<=0 ) RETURN
  kk = Nz + 1
  s1 = Y(kk)
  Y(kk) = s1*csr(1)
  IF ( inu==1 ) RETURN
  kk = Nz + 2
  s2 = Y(kk)
  Y(kk) = s2*csr(1)
  IF ( inu==2 ) RETURN
  t2 = Fnu + (kk-1)
  ck = CMPLX(t2,0.0E0)*rz
  kflag = 1
  GOTO 500
  !-----------------------------------------------------------------------
  !     FNU=HALF ODD INTEGER CASE, DNU=-0.5
  !-----------------------------------------------------------------------
  800  s1 = coef
  s2 = coef
  GOTO 200
  99999 CONTINUE
  END SUBROUTINE CBKNU
