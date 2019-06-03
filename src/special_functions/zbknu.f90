!** ZBKNU
SUBROUTINE ZBKNU(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Tol,Elim,Alim)
  !>
  !  Subsidiary to ZAIRY, ZBESH, ZBESI and ZBESK
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
  USE service, ONLY : D1MACH, I1MACH
  !
  INTEGER i, iflag, inu, k, kflag, kk, Kode, koded, N, Nz, idum, j, ic, inub, nw
  REAL(DP) :: aa, ak, Alim, ascle, a1, a2, bb, bk, bry(3), caz, cbi, cbr, cchi, &
    cchr, cki, ckr, coefi, coefr, crscr, csclr, cshi, cshr, csi, csr, csrr(3), &
    cssr(3), czi, czr, dnu, dnu2, Elim, etest, fc, fhs, fi, fk, fks, fmui, fmur, &
    Fnu, fr, g1, g2, pi, pr, pti, ptr, p1i, p1r, p2i, p2m, p2r, qi, qr, rak, &
    rcaz, rzi, rzr, s, smui, smur, sti, str, s1i, s1r, s2i, s2r, tm, Tol, t1, &
    t2, Yi(N), Yr(N), Zi, Zr, elm, celmr, zdr, zdi, as, alas, helim, cyr(2), cyi(2)
  !     COMPLEX Z,Y,A,B,RZ,SMU,FU,FMU,F,FLRZ,CZ,S1,S2,CSH,CCH
  !     COMPLEX CK,P,Q,COEF,P1,P2,CBK,PT,CZERO,CONE,CTWO,ST,EZ,CS,DK
  !
  INTEGER, PARAMETER :: kmax = 30
  REAL(DP), PARAMETER :: czeror = 0.0D0, czeroi = 0.0D0, coner = 1.0D0, &
    conei = 0.0D0, ctwor = 2.0D0, r1 = 2.0D0
  REAL(DP), PARAMETER :: dpi = 3.14159265358979324D0, rthpi= 1.25331413731550025D0, &
    spi = 1.90985931710274403D0, hpi = 1.57079632679489662D0, &
    fpi = 1.89769999331517738D0, tth = 6.66666666666666666D-01
  REAL(DP), PARAMETER :: cc(8) = [ 5.77215664901532861D-01, -4.20026350340952355D-02, &
    -4.21977345555443367D-02, 7.21894324666309954D-03,-2.15241674114950973D-04, &
    -2.01348547807882387D-05, 1.13302723198169588D-06, 6.11609510448141582D-09 ]
  !* FIRST EXECUTABLE STATEMENT  ZBKNU
  caz = ZABS(Zr,Zi)
  csclr = 1.0D0/Tol
  crscr = Tol
  cssr(1) = csclr
  cssr(2) = 1.0D0
  cssr(3) = crscr
  csrr(1) = crscr
  csrr(2) = 1.0D0
  csrr(3) = csclr
  bry(1) = 1.0D+3*D1MACH(1)/Tol
  bry(2) = 1.0D0/bry(1)
  bry(3) = D1MACH(2)
  Nz = 0
  iflag = 0
  koded = Kode
  rcaz = 1.0D0/caz
  str = Zr*rcaz
  sti = -Zi*rcaz
  rzr = (str+str)*rcaz
  rzi = (sti+sti)*rcaz
  inu = INT( Fnu + 0.5D0 )
  dnu = Fnu - inu
  IF ( ABS(dnu)/=0.5D0 ) THEN
    dnu2 = 0.0D0
    IF ( ABS(dnu)>Tol ) dnu2 = dnu*dnu
    IF ( caz<=r1 ) THEN
      !-----------------------------------------------------------------------
      !     SERIES FOR ABS(Z).LE.R1
      !-----------------------------------------------------------------------
      fc = 1.0D0
      CALL ZLOG(rzr,rzi,smur,smui,idum)
      fmur = smur*dnu
      fmui = smui*dnu
      CALL ZSHCH(fmur,fmui,cshr,cshi,cchr,cchi)
      IF ( dnu/=0.0D0 ) THEN
        fc = dnu*dpi
        fc = fc/SIN(fc)
        smur = cshr/dnu
        smui = cshi/dnu
      END IF
      a2 = 1.0D0 + dnu
      !-----------------------------------------------------------------------
      !     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
      !-----------------------------------------------------------------------
      t2 = EXP(-DGAMLN(a2,idum))
      t1 = 1.0D0/(t2*fc)
      IF ( ABS(dnu)>0.1D0 ) THEN
        g1 = (t1-t2)/(dnu+dnu)
      ELSE
        !-----------------------------------------------------------------------
        !     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
        !-----------------------------------------------------------------------
        ak = 1.0D0
        s = cc(1)
        DO k = 2, 8
          ak = ak*dnu2
          tm = cc(k)*ak
          s = s + tm
          IF ( ABS(tm)<Tol ) EXIT
        END DO
        g1 = -s
      END IF
      g2 = (t1+t2)*0.5D0
      fr = fc*(cchr*g1+smur*g2)
      fi = fc*(cchi*g1+smui*g2)
      CALL ZEXP(fmur,fmui,str,sti)
      pr = 0.5D0*str/t2
      pi = 0.5D0*sti/t2
      CALL ZDIV(0.5D0,0.0D0,str,sti,ptr,pti)
      qr = ptr/t1
      qi = pti/t1
      s1r = fr
      s1i = fi
      s2r = pr
      s2i = pi
      ak = 1.0D0
      a1 = 1.0D0
      ckr = coner
      cki = conei
      bk = 1.0D0 - dnu2
      IF ( inu>0.OR.N>1 ) THEN
        !-----------------------------------------------------------------------
        !     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
        !-----------------------------------------------------------------------
        IF ( caz>=Tol ) THEN
          CALL ZMLT(Zr,Zi,Zr,Zi,czr,czi)
          czr = 0.25D0*czr
          czi = 0.25D0*czi
          t1 = 0.25D0*caz*caz
          DO
            fr = (fr*ak+pr+qr)/bk
            fi = (fi*ak+pi+qi)/bk
            str = 1.0D0/(ak-dnu)
            pr = pr*str
            pi = pi*str
            str = 1.0D0/(ak+dnu)
            qr = qr*str
            qi = qi*str
            str = ckr*czr - cki*czi
            rak = 1.0D0/ak
            cki = (ckr*czi+cki*czr)*rak
            ckr = str*rak
            s1r = ckr*fr - cki*fi + s1r
            s1i = ckr*fi + cki*fr + s1i
            str = pr - fr*ak
            sti = pi - fi*ak
            s2r = ckr*str - cki*sti + s2r
            s2i = ckr*sti + cki*str + s2i
            a1 = a1*t1*rak
            bk = bk + ak + ak + 1.0D0
            ak = ak + 1.0D0
            IF ( a1<=Tol ) EXIT
          END DO
        END IF
        kflag = 2
        a1 = Fnu + 1.0D0
        ak = a1*ABS(smur)
        IF ( ak>Alim ) kflag = 3
        str = cssr(kflag)
        p2r = s2r*str
        p2i = s2i*str
        CALL ZMLT(p2r,p2i,rzr,rzi,s2r,s2i)
        s1r = s1r*str
        s1i = s1i*str
        IF ( koded/=1 ) THEN
          CALL ZEXP(Zr,Zi,fr,fi)
          CALL ZMLT(s1r,s1i,fr,fi,s1r,s1i)
          CALL ZMLT(s2r,s2i,fr,fi,s2r,s2i)
        END IF
        GOTO 200
      ELSE
        !-----------------------------------------------------------------------
        !     GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
        !-----------------------------------------------------------------------
        IF ( caz>=Tol ) THEN
          CALL ZMLT(Zr,Zi,Zr,Zi,czr,czi)
          czr = 0.25D0*czr
          czi = 0.25D0*czi
          t1 = 0.25D0*caz*caz
          DO
            fr = (fr*ak+pr+qr)/bk
            fi = (fi*ak+pi+qi)/bk
            str = 1.0D0/(ak-dnu)
            pr = pr*str
            pi = pi*str
            str = 1.0D0/(ak+dnu)
            qr = qr*str
            qi = qi*str
            str = ckr*czr - cki*czi
            rak = 1.0D0/ak
            cki = (ckr*czi+cki*czr)*rak
            ckr = str*rak
            s1r = ckr*fr - cki*fi + s1r
            s1i = ckr*fi + cki*fr + s1i
            a1 = a1*t1*rak
            bk = bk + ak + ak + 1.0D0
            ak = ak + 1.0D0
            IF ( a1<=Tol ) EXIT
          END DO
        END IF
        Yr(1) = s1r
        Yi(1) = s1i
        IF ( koded==1 ) RETURN
        CALL ZEXP(Zr,Zi,str,sti)
        CALL ZMLT(s1r,s1i,str,sti,Yr(1),Yi(1))
        RETURN
      END IF
    END IF
  END IF
  !-----------------------------------------------------------------------
  !     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
  !     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
  !     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
  !     RECURSION
  !-----------------------------------------------------------------------
  CALL ZSQRT(Zr,Zi,str,sti)
  CALL ZDIV(rthpi,czeroi,str,sti,coefr,coefi)
  kflag = 2
  IF ( koded/=2 ) THEN
    IF ( Zr>Alim ) THEN
      !-----------------------------------------------------------------------
      !     SCALE BY EXP(Z), IFLAG = 1 CASES
      !-----------------------------------------------------------------------
      koded = 2
      iflag = 1
      kflag = 2
    ELSE
      !     BLANK LINE
      str = EXP(-Zr)*cssr(kflag)
      sti = -str*SIN(Zi)
      str = str*COS(Zi)
      CALL ZMLT(coefr,coefi,str,sti,coefr,coefi)
    END IF
  END IF
  IF ( ABS(dnu)==0.5D0 ) GOTO 800
  !-----------------------------------------------------------------------
  !     MILLER ALGORITHM FOR ABS(Z).GT.R1
  !-----------------------------------------------------------------------
  ak = COS(dpi*dnu)
  ak = ABS(ak)
  IF ( ak==czeror ) GOTO 800
  fhs = ABS(0.25D0-dnu2)
  IF ( fhs==czeror ) GOTO 800
  !-----------------------------------------------------------------------
  !     COMPUTE R2=F(E). IF ABS(Z).GE.R2, USE FORWARD RECURRENCE TO
  !     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
  !     12.LE.E.LE.60. E IS COMPUTED FROM 2**(-E)=B**(1-I1MACH(14))=
  !     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
  !-----------------------------------------------------------------------
  t1 = I1MACH(14) - 1
  t1 = t1*D1MACH(5)*3.321928094D0
  t1 = MAX(t1,12.0D0)
  t1 = MIN(t1,60.0D0)
  t2 = tth*t1 - 6.0D0
  IF ( Zr/=0.0D0 ) THEN
    t1 = ATAN(Zi/Zr)
    t1 = ABS(t1)
  ELSE
    t1 = hpi
  END IF
  IF ( t2>caz ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE BACKWARD INDEX K FOR ABS(Z).LT.R2
    !-----------------------------------------------------------------------
    a2 = SQRT(caz)
    ak = fpi*ak/(Tol*SQRT(a2))
    aa = 3.0D0*t1/(1.0D0+caz)
    bb = 14.7D0*t1/(28.0D0+caz)
    ak = (LOG(ak)+caz*COS(aa)/(1.0D0+0.008D0*caz))/COS(bb)
    fk = 0.12125D0*ak*ak/caz + 1.5D0
  ELSE
    !-----------------------------------------------------------------------
    !     FORWARD RECURRENCE LOOP WHEN ABS(Z).GE.R2
    !-----------------------------------------------------------------------
    etest = ak/(dpi*caz*Tol)
    fk = coner
    IF ( etest<coner ) GOTO 100
    fks = ctwor
    ckr = caz + caz + ctwor
    p1r = czeror
    p2r = coner
    DO i = 1, kmax
      ak = fhs/fks
      cbr = ckr/(fk+coner)
      ptr = p2r
      p2r = cbr*p2r - p1r*ak
      p1r = ptr
      ckr = ckr + ctwor
      fks = fks + fk + fk + ctwor
      fhs = fhs + fk + fk
      fk = fk + coner
      str = ABS(p2r)*fk
      IF ( etest<str ) GOTO 50
    END DO
    !
    !
    Nz = -2
    RETURN
    50  fk = fk + spi*t1*SQRT(t2/caz)
    fhs = ABS(0.25D0-dnu2)
  END IF
  !-----------------------------------------------------------------------
  !     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
  !-----------------------------------------------------------------------
  100  k = INT( fk )
  fk = k
  fks = fk*fk
  p1r = czeror
  p1i = czeroi
  p2r = Tol
  p2i = czeroi
  csr = p2r
  csi = p2i
  DO i = 1, k
    a1 = fks - fk
    ak = (fks+fk)/(a1+fhs)
    rak = 2.0D0/(fk+coner)
    cbr = (fk+Zr)*rak
    cbi = Zi*rak
    ptr = p2r
    pti = p2i
    p2r = (ptr*cbr-pti*cbi-p1r)*ak
    p2i = (pti*cbr+ptr*cbi-p1i)*ak
    p1r = ptr
    p1i = pti
    csr = csr + p2r
    csi = csi + p2i
    fks = a1 - fk + coner
    fk = fk - coner
  END DO
  !-----------------------------------------------------------------------
  !     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER
  !     SCALING
  !-----------------------------------------------------------------------
  tm = ZABS(csr,csi)
  ptr = 1.0D0/tm
  s1r = p2r*ptr
  s1i = p2i*ptr
  csr = csr*ptr
  csi = -csi*ptr
  CALL ZMLT(coefr,coefi,s1r,s1i,str,sti)
  CALL ZMLT(str,sti,csr,csi,s1r,s1i)
  IF ( inu>0.OR.N>1 ) THEN
    !-----------------------------------------------------------------------
    !     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
    !-----------------------------------------------------------------------
    tm = ZABS(p2r,p2i)
    ptr = 1.0D0/tm
    p1r = p1r*ptr
    p1i = p1i*ptr
    p2r = p2r*ptr
    p2i = -p2i*ptr
    CALL ZMLT(p1r,p1i,p2r,p2i,ptr,pti)
    str = dnu + 0.5D0 - ptr
    sti = -pti
    CALL ZDIV(str,sti,Zr,Zi,str,sti)
    str = str + 1.0D0
    CALL ZMLT(str,sti,s1r,s1i,s2r,s2i)
  ELSE
    zdr = Zr
    zdi = Zi
    IF ( iflag/=1 ) GOTO 400
    GOTO 700
  END IF
  !-----------------------------------------------------------------------
  !     FORWARD RECURSION ON THE THREE TERM RECURSION WITH RELATION WITH
  !     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
  !-----------------------------------------------------------------------
  200  str = dnu + 1.0D0
  ckr = str*rzr
  cki = str*rzi
  IF ( N==1 ) inu = inu - 1
  IF ( inu>0 ) THEN
    inub = 1
    IF ( iflag==1 ) THEN
      !-----------------------------------------------------------------------
      !     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
      !-----------------------------------------------------------------------
      helim = 0.5D0*Elim
      elm = EXP(-Elim)
      celmr = elm
      ascle = bry(1)
      zdr = Zr
      zdi = Zi
      ic = -1
      j = 2
      DO i = 1, inu
        str = s2r
        sti = s2i
        s2r = str*ckr - sti*cki + s1r
        s2i = sti*ckr + str*cki + s1i
        s1r = str
        s1i = sti
        ckr = ckr + rzr
        cki = cki + rzi
        as = ZABS(s2r,s2i)
        alas = LOG(as)
        p2r = -zdr + alas
        IF ( p2r>=(-Elim) ) THEN
          CALL ZLOG(s2r,s2i,str,sti,idum)
          p2r = -zdr + str
          p2i = -zdi + sti
          p2m = EXP(p2r)/Tol
          p1r = p2m*COS(p2i)
          p1i = p2m*SIN(p2i)
          CALL ZUCHK(p1r,p1i,nw,ascle,Tol)
          IF ( nw==0 ) THEN
            j = 3 - j
            cyr(j) = p1r
            cyi(j) = p1i
            IF ( ic==(i-1) ) GOTO 600
            ic = i
            CYCLE
          END IF
        END IF
        IF ( alas>=helim ) THEN
          zdr = zdr - Elim
          s1r = s1r*celmr
          s1i = s1i*celmr
          s2r = s2r*celmr
          s2i = s2i*celmr
        END IF
      END DO
      IF ( N==1 ) THEN
        s1r = s2r
        s1i = s2i
      END IF
      GOTO 700
    END IF
  ELSE
    IF ( N<=1 ) THEN
      s1r = s2r
      s1i = s2i
    END IF
    zdr = Zr
    zdi = Zi
    IF ( iflag/=1 ) GOTO 400
    GOTO 700
  END IF
  300  p1r = csrr(kflag)
  ascle = bry(kflag)
  DO i = inub, inu
    str = s2r
    sti = s2i
    s2r = ckr*str - cki*sti + s1r
    s2i = ckr*sti + cki*str + s1i
    s1r = str
    s1i = sti
    ckr = ckr + rzr
    cki = cki + rzi
    IF ( kflag<3 ) THEN
      p2r = s2r*p1r
      p2i = s2i*p1r
      str = ABS(p2r)
      sti = ABS(p2i)
      p2m = MAX(str,sti)
      IF ( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*p1r
        s1i = s1i*p1r
        s2r = p2r
        s2i = p2i
        str = cssr(kflag)
        s1r = s1r*str
        s1i = s1i*str
        s2r = s2r*str
        s2i = s2i*str
        p1r = csrr(kflag)
      END IF
    END IF
  END DO
  IF ( N==1 ) THEN
    s1r = s2r
    s1i = s2i
  END IF
  400  str = csrr(kflag)
  Yr(1) = s1r*str
  Yi(1) = s1i*str
  IF ( N==1 ) RETURN
  Yr(2) = s2r*str
  Yi(2) = s2i*str
  IF ( N==2 ) RETURN
  kk = 2
  500  kk = kk + 1
  IF ( kk>N ) RETURN
  p1r = csrr(kflag)
  ascle = bry(kflag)
  DO i = kk, N
    p2r = s2r
    p2i = s2i
    s2r = ckr*p2r - cki*p2i + s1r
    s2i = cki*p2r + ckr*p2i + s1i
    s1r = p2r
    s1i = p2i
    ckr = ckr + rzr
    cki = cki + rzi
    p2r = s2r*p1r
    p2i = s2i*p1r
    Yr(i) = p2r
    Yi(i) = p2i
    IF ( kflag<3 ) THEN
      str = ABS(p2r)
      sti = ABS(p2i)
      p2m = MAX(str,sti)
      IF ( p2m>ascle ) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1r = s1r*p1r
        s1i = s1i*p1r
        s2r = p2r
        s2i = p2i
        str = cssr(kflag)
        s1r = s1r*str
        s1i = s1i*str
        s2r = s2r*str
        s2i = s2i*str
        p1r = csrr(kflag)
      END IF
    END IF
  END DO
  RETURN
  600  kflag = 1
  inub = i + 1
  s2r = cyr(j)
  s2i = cyi(j)
  j = 3 - j
  s1r = cyr(j)
  s1i = cyi(j)
  IF ( inub<=inu ) GOTO 300
  IF ( N==1 ) THEN
    s1r = s2r
    s1i = s2i
  END IF
  GOTO 400
  700  Yr(1) = s1r
  Yi(1) = s1i
  IF ( N/=1 ) THEN
    Yr(2) = s2r
    Yi(2) = s2i
  END IF
  ascle = bry(1)
  CALL ZKSCL(zdr,zdi,Fnu,N,Yr,Yi,Nz,rzr,rzi,ascle,Tol,Elim)
  inu = N - Nz
  IF ( inu<=0 ) RETURN
  kk = Nz + 1
  s1r = Yr(kk)
  s1i = Yi(kk)
  Yr(kk) = s1r*csrr(1)
  Yi(kk) = s1i*csrr(1)
  IF ( inu==1 ) RETURN
  kk = Nz + 2
  s2r = Yr(kk)
  s2i = Yi(kk)
  Yr(kk) = s2r*csrr(1)
  Yi(kk) = s2i*csrr(1)
  IF ( inu==2 ) RETURN
  t2 = Fnu + (kk-1)
  ckr = t2*rzr
  cki = t2*rzi
  kflag = 1
  GOTO 500
  !-----------------------------------------------------------------------
  !     FNU=HALF ODD INTEGER CASE, DNU=-0.5
  !-----------------------------------------------------------------------
  800  s1r = coefr
  s1i = coefi
  s2r = coefr
  s2i = coefi
  GOTO 200
  RETURN
END SUBROUTINE ZBKNU
