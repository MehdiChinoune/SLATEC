!** ZBKNU
PURE SUBROUTINE ZBKNU(Z,Fnu,Kode,N,Y,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZAIRY, ZBESH, ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBKNU-A, ZBKNU-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE.
  !
  !***
  ! **See also:**  ZAIRY, ZBESH, ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, DGAMLN, I1MACH, ZABS, ZDIV, ZEXP, ZKSCL,
  !                    ZLOG, ZMLT, ZSHCH, ZSQRT, ZUCHK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP, ZLOG and ZSQRT to EXTERNAL statement.  (RWC)
  USE service, ONLY : tiny_dp, huge_dp, log10_radix_dp, digits_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, iflag, inu, k, kflag, kk, koded, nw, j, ic, inub
  COMPLEX(DP) :: cch, ck, coef, crsc, cs, cscl, csh, csr(3), css(3), cz, f, fmu, p, &
    pt, p1, p2, q, rz, smu, st, s1, s2, zd, celm, cy(2)
  REAL(DP) :: aa, ak, ascle, a1, a2, bb, bk, bry(3), caz, dnu, dnu2, etest, fc, &
    fhs, fk, fks, g1, g2, p2i, p2m, p2r, rk, s, tm, t1, t2, xx, yy, helim, elm, &
    xd, yd, alas, as
  !
  INTEGER, PARAMETER :: kmax = 30
  REAL(DP), PARAMETER ::  r1 = 2._DP
  !
  REAL(DP), PARAMETER ::  pi = 3.14159265358979324_DP, rthpi = 1.25331413731550025_DP, &
    spi = 1.90985931710274403_DP, hpi = 1.57079632679489662_DP, &
    fpi = 1.89769999331517738_DP, tth = 6.66666666666666666E-01_DP
  !
  REAL(DP), PARAMETER :: cc(8) = [ 5.77215664901532861E-01_DP, -4.20026350340952355E-02_DP, &
    -4.21977345555443367E-02_DP, 7.21894324666309954E-03_DP, -2.15241674114950973E-04_DP, &
    -2.01348547807882387E-05_DP, 1.13302723198169588E-06_DP, 6.11609510448141582E-09_DP ]
  !
  !* FIRST EXECUTABLE STATEMENT  ZBKNU
  xx = REAL(Z,DP)
  yy = AIMAG(Z)
  caz = ABS(Z)
  cscl = CMPLX(1._DP/Tol,0._DP,DP)
  crsc = CMPLX(Tol,0._DP,DP)
  css(1) = cscl
  css(2) = (1._DP,0._DP)
  css(3) = crsc
  csr(1) = crsc
  csr(2) = (1._DP,0._DP)
  csr(3) = cscl
  bry(1) = 1.E+3_DP*tiny_dp/Tol
  bry(2) = 1._DP/bry(1)
  bry(3) = huge_dp
  Nz = 0
  iflag = 0
  koded = Kode
  rz = (2._DP,0._DP)/Z
  inu = INT( Fnu + 0.5_DP )
  dnu = Fnu - inu
  IF( ABS(dnu)/=0.5_DP ) THEN
    dnu2 = 0._DP
    IF( ABS(dnu)>Tol ) dnu2 = dnu*dnu
    IF( caz<=r1 ) THEN
      !-----------------------------------------------------------------------
      !     SERIES FOR ABS(Z)<=R1
      !-----------------------------------------------------------------------
      fc = 1._DP
      smu = LOG(rz)
      fmu = smu*CMPLX(dnu,0._DP,DP)
      CALL ZSHCH(fmu,csh,cch)
      IF( dnu/=0._DP ) THEN
        fc = dnu*pi
        fc = fc/SIN(fc)
        smu = csh*CMPLX(1._DP/dnu,0._DP,DP)
      END IF
      a2 = 1._DP + dnu
      !-----------------------------------------------------------------------
      !     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
      !-----------------------------------------------------------------------
      t2 = EXP(-LOG_GAMMA(a2))
      t1 = 1._DP/(t2*fc)
      IF( ABS(dnu)>0.1_DP ) THEN
        g1 = (t1-t2)/(dnu+dnu)
      ELSE
        !-----------------------------------------------------------------------
        !     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
        !-----------------------------------------------------------------------
        ak = 1._DP
        s = cc(1)
        DO k = 2, 8
          ak = ak*dnu2
          tm = cc(k)*ak
          s = s + tm
          IF( ABS(tm)<Tol ) EXIT
        END DO
        g1 = -s
      END IF
      g2 = 0.5_DP*(t1+t2)*fc
      g1 = g1*fc
      f = CMPLX(g1,0._DP,DP)*cch + smu*CMPLX(g2,0._DP,DP)
      pt = EXP(fmu)
      p = CMPLX(0.5_DP/t2,0._DP,DP)*pt
      q = CMPLX(0.5_DP/t1,0._DP,DP)/pt
      s1 = f
      s2 = p
      ak = 1._DP
      a1 = 1._DP
      ck = (1._DP,0._DP)
      bk = 1._DP - dnu2
      IF( inu>0 .OR. N>1 ) THEN
        !-----------------------------------------------------------------------
        !     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
        !-----------------------------------------------------------------------
        IF( caz>=Tol ) THEN
          cz = Z*Z*CMPLX(0.25_DP,0._DP,DP)
          t1 = 0.25_DP*caz*caz
          DO
            f = (f*CMPLX(ak,0._DP,DP)+p+q)*CMPLX(1._DP/bk,0._DP,DP)
            p = p*CMPLX(1._DP/(ak-dnu),0._DP,DP)
            q = q*CMPLX(1._DP/(ak+dnu),0._DP,DP)
            rk = 1._DP/ak
            ck = ck*cz*CMPLX(rk,0._DP,DP)
            s1 = s1 + ck*f
            s2 = s2 + ck*(p-f*CMPLX(ak,0._DP,DP))
            a1 = a1*t1*rk
            bk = bk + ak + ak + 1._DP
            ak = ak + 1._DP
            IF( a1<=Tol ) EXIT
          END DO
        END IF
        kflag = 2
        bk = REAL(smu,DP)
        a1 = Fnu + 1._DP
        ak = a1*ABS(bk)
        IF( ak>Alim ) kflag = 3
        p2 = s2*css(kflag)
        s2 = p2*rz
        s1 = s1*css(kflag)
        IF( koded/=1 ) THEN
          f = EXP(Z)
          s1 = s1*f
          s2 = s2*f
        END IF
        GOTO 200
      ELSE
        !-----------------------------------------------------------------------
        !     GENERATE K(FNU,Z), 0.0D0 <= FNU < 0.5D0 AND N=1
        !-----------------------------------------------------------------------
        IF( caz>=Tol ) THEN
          cz = Z*Z*CMPLX(0.25_DP,0._DP,DP)
          t1 = 0.25_DP*caz*caz
          DO
            f = (f*CMPLX(ak,0._DP,DP)+p+q)*CMPLX(1._DP/bk,0._DP,DP)
            p = p*CMPLX(1._DP/(ak-dnu),0._DP,DP)
            q = q*CMPLX(1._DP/(ak+dnu),0._DP,DP)
            rk = 1._DP/ak
            ck = ck*cz*CMPLX(rk,0._DP,DP)
            s1 = s1 + ck*f
            a1 = a1*t1*rk
            bk = bk + ak + ak + 1._DP
            ak = ak + 1._DP
            IF( a1<=Tol ) EXIT
          END DO
        END IF
        Y(1) = s1
        IF( koded==1 ) RETURN
        Y(1) = s1*EXP(Z)
        RETURN
      END IF
    END IF
  END IF
  !-----------------------------------------------------------------------
  !     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
  !     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
  !     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD RECURSION
  !-----------------------------------------------------------------------
  coef = CMPLX(rthpi,0._DP,DP)/SQRT(Z)
  kflag = 2
  IF( koded/=2 ) THEN
    IF( xx>Alim ) THEN
      !-----------------------------------------------------------------------
      !     SCALE BY EXP(Z), IFLAG = 1 CASES
      !-----------------------------------------------------------------------
      koded = 2
      iflag = 1
      kflag = 2
    ELSE
      !     BLANK LINE
      a1 = EXP(-xx)*REAL(css(kflag),DP)
      pt = CMPLX(a1,0._DP,DP)*CMPLX(COS(yy),-SIN(yy),DP)
      coef = coef*pt
    END IF
  END IF
  IF( ABS(dnu)==0.5_DP ) GOTO 800
  !-----------------------------------------------------------------------
  !     MILLER ALGORITHM FOR ABS(Z)>R1
  !-----------------------------------------------------------------------
  ak = COS(pi*dnu)
  ak = ABS(ak)
  IF( ak==0._DP ) GOTO 800
  fhs = ABS(0.25_DP-dnu2)
  IF( fhs==0._DP ) GOTO 800
  !-----------------------------------------------------------------------
  !     COMPUTE R2=F(E). IF ABS(Z)>=R2, USE FORWARD RECURRENCE TO
  !     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
  !     12<=E<=60. E IS COMPUTED FROM 2**(-E)=B**(1-digits_dp)=
  !     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
  !-----------------------------------------------------------------------
  t1 = (digits_dp-1)*log10_radix_dp*3.321928094_DP
  t1 = MAX(t1,12._DP)
  t1 = MIN(t1,60._DP)
  t2 = tth*t1 - 6._DP
  IF( xx/=0._DP ) THEN
    t1 = ATAN(yy/xx)
    t1 = ABS(t1)
  ELSE
    t1 = hpi
  END IF
  IF( t2>caz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE BACKWARD INDEX K FOR ABS(Z)<R2
    !-----------------------------------------------------------------------
    a2 = SQRT(caz)
    ak = fpi*ak/(Tol*SQRT(a2))
    aa = 3._DP*t1/(1._DP+caz)
    bb = 14.7_DP*t1/(28._DP+caz)
    ak = (LOG(ak)+caz*COS(aa)/(1._DP+0.008_DP*caz))/COS(bb)
    fk = 0.12125_DP*ak*ak/caz + 1.5_DP
  ELSE
    !-----------------------------------------------------------------------
    !     FORWARD RECURRENCE LOOP WHEN ABS(Z)>=R2
    !-----------------------------------------------------------------------
    etest = ak/(pi*caz*Tol)
    fk = 1._DP
    IF( etest<1._DP ) GOTO 100
    fks = 2._DP
    rk = caz + caz + 2._DP
    a1 = 0._DP
    a2 = 1._DP
    DO i = 1, kmax
      ak = fhs/fks
      bk = rk/(fk+1._DP)
      tm = a2
      a2 = bk*a2 - ak*a1
      a1 = tm
      rk = rk + 2._DP
      fks = fks + fk + fk + 2._DP
      fhs = fhs + fk + fk
      fk = fk + 1._DP
      tm = ABS(a2)*fk
      IF( etest<tm ) GOTO 50
    END DO
    Nz = -2
    RETURN
    50  fk = fk + spi*t1*SQRT(t2/caz)
    fhs = ABS(0.25_DP-dnu2)
  END IF
  100  k = INT( fk )
  !-----------------------------------------------------------------------
  !     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  fk = k
  fks = fk*fk
  p1 = (0._DP,0._DP)
  p2 = CMPLX(Tol,0._DP,DP)
  cs = p2
  DO i = 1, k
    a1 = fks - fk
    a2 = (fks+fk)/(a1+fhs)
    rk = 2._DP/(fk+1._DP)
    t1 = (fk+xx)*rk
    t2 = yy*rk
    pt = p2
    p2 = (p2*CMPLX(t1,t2,DP)-p1)*CMPLX(a2,0._DP,DP)
    p1 = pt
    cs = cs + p2
    fks = a1 - fk + 1._DP
    fk = fk - 1._DP
  END DO
  !-----------------------------------------------------------------------
  !     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER SCALING
  !-----------------------------------------------------------------------
  tm = ABS(cs)
  pt = CMPLX(1._DP/tm,0._DP,DP)
  s1 = pt*p2
  cs = CONJG(cs)*pt
  s1 = coef*s1*cs
  IF( inu>0 .OR. N>1 ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
    !-----------------------------------------------------------------------
    tm = ABS(p2)
    pt = CMPLX(1._DP/tm,0._DP,DP)
    p1 = pt*p1
    p2 = CONJG(p2)*pt
    pt = p1*p2
    s2 = s1*((1._DP,0._DP)+(CMPLX(dnu+0.5_DP,0._DP,DP)-pt)/Z)
  ELSE
    zd = Z
    IF( iflag/=1 ) GOTO 400
    GOTO 700
  END IF
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
  !     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
  !-----------------------------------------------------------------------
  200  ck = CMPLX(dnu+1._DP,0._DP,DP)*rz
  IF( N==1 ) inu = inu - 1
  IF( inu>0 ) THEN
    inub = 1
    IF( iflag==1 ) THEN
      !-----------------------------------------------------------------------
      !     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
      !-----------------------------------------------------------------------
      helim = 0.5_DP*Elim
      elm = EXP(-Elim)
      celm = CMPLX(elm,0._DP,DP)
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
        alas = LOG(as)
        p2r = -xd + alas
        IF( p2r>=(-Elim) ) THEN
          p2 = -zd + LOG(s2)
          p2r = REAL(p2,DP)
          p2i = AIMAG(p2)
          p2m = EXP(p2r)/Tol
          p1 = CMPLX(p2m,0._DP,DP)*CMPLX(COS(p2i),SIN(p2i),DP)
          CALL ZUCHK(p1,nw,ascle,Tol)
          IF( nw==0 ) THEN
            j = 3 - j
            cy(j) = p1
            IF( ic==(i-1) ) GOTO 600
            ic = i
            CYCLE
          END IF
        END IF
        IF( alas>=helim ) THEN
          xd = xd - Elim
          s1 = s1*celm
          s2 = s2*celm
          zd = CMPLX(xd,yd,DP)
        END IF
      END DO
      IF( N==1 ) s1 = s2
      GOTO 700
    END IF
  ELSE
    IF( N==1 ) s1 = s2
    zd = Z
    IF( iflag/=1 ) GOTO 400
    GOTO 700
  END IF
  300  p1 = csr(kflag)
  ascle = bry(kflag)
  DO i = inub, inu
    st = s2
    s2 = ck*s2 + s1
    s1 = st
    ck = ck + rz
    IF( kflag<3 ) THEN
      p2 = s2*p1
      p2r = REAL(p2,DP)
      p2i = AIMAG(p2)
      p2r = ABS(p2r)
      p2i = ABS(p2i)
      p2m = MAX(p2r,p2i)
      IF( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
      END IF
    END IF
  END DO
  IF( N==1 ) s1 = s2
  400  Y(1) = s1*csr(kflag)
  IF( N==1 ) RETURN
  Y(2) = s2*csr(kflag)
  IF( N==2 ) RETURN
  kk = 2
  500  kk = kk + 1
  IF( kk>N ) RETURN
  p1 = csr(kflag)
  ascle = bry(kflag)
  DO i = kk, N
    p2 = s2
    s2 = ck*s2 + s1
    s1 = p2
    ck = ck + rz
    p2 = s2*p1
    Y(i) = p2
    IF( kflag<3 ) THEN
      p2r = REAL(p2,DP)
      p2i = AIMAG(p2)
      p2r = ABS(p2r)
      p2i = ABS(p2i)
      p2m = MAX(p2r,p2i)
      IF( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1*p1
        s2 = p2
        s1 = s1*css(kflag)
        s2 = s2*css(kflag)
        p1 = csr(kflag)
      END IF
    END IF
  END DO
  RETURN
  600  kflag = 1
  inub = i + 1
  s2 = cy(j)
  j = 3 - j
  s1 = cy(j)
  IF( inub<=inu ) GOTO 300
  IF( N==1 ) s1 = s2
  GOTO 400
  700  Y(1) = s1
  IF( N/=1 ) Y(2) = s2
  ascle = bry(1)
  CALL ZKSCL(zd,Fnu,N,Y,Nz,rz,ascle,Tol,Elim)
  inu = N - Nz
  IF( inu<=0 ) RETURN
  kk = Nz + 1
  s1 = Y(kk)
  Y(kk) = s1*csr(1)
  IF( inu==1 ) RETURN
  kk = Nz + 2
  s2 = Y(kk)
  Y(kk) = s2*csr(1)
  IF( inu==2 ) RETURN
  t2 = Fnu + (kk-1)
  ck = CMPLX(t2,0._DP,DP)*rz
  kflag = 1
  GOTO 500
  !-----------------------------------------------------------------------
  !     FNU=HALF ODD INTEGER CASE, DNU=-0.5
  !-----------------------------------------------------------------------
  800  s1 = coef
  s2 = coef
  GOTO 200
  !
  RETURN
END SUBROUTINE ZBKNU