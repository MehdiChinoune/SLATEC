!*==ZUNI2.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZUNI2
SUBROUTINE ZUNI2(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Nlast,Fnul,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--ZUNI25
  !***BEGIN PROLOGUE  ZUNI2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUNI2-A, ZUNI2-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
  !     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
  !     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
  !
  !     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
  !     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
  !     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
  !     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
  !     Y(I)=CZERO FOR I=NLAST+1,N
  !
  !***SEE ALSO  ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, ZABS, ZAIRY, ZUCHK, ZUNHJ, ZUOIK
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZUNI2
  !     COMPLEX AI,ARG,ASUM,BSUM,CFN,CI,CID,CIP,CONE,CRSC,CSCL,CSR,CSS,
  !    *CZERO,C1,C2,DAI,PHI,RZ,S1,S2,Y,Z,ZB,ZETA1,ZETA2,ZN
  REAL(8) :: aarg, aic, aii, air, Alim, ang, aphi, argi, &
    argr, ascle, asumi, asumr, bry, bsumi, bsumr, &
    cidi, cipi, cipr, coner, crsc, cscl, csrr, cssr, &
    c1r, c2i, c2m, c2r, daii, dair, Elim, fn, Fnu, &
    Fnul, hpi, phii, phir, rast, raz, rs1, rzi, rzr, &
    sti, str, s1i, s1r, s2i, s2r, Tol, Yi, Yr, zbi, &
    zbr, zeroi, zeror, zeta1i, zeta1r, zeta2i, zeta2r, &
    Zi, zni, znr, Zr, cyr, cyi, D1MACH, ZABS, car, &
    sar
  INTEGER i, iflag, in, inu, j, k, Kode, N, nai, nd, ndai, &
    Nlast, nn, nuf, nw, Nz, idum
  DIMENSION bry(3), Yr(N), Yi(N), cipr(4), cipi(4), cssr(3), csrr(3), &
    cyr(2), cyi(2)
  EXTERNAL ZABS
  DATA zeror, zeroi, coner/0.0D0, 0.0D0, 1.0D0/
  DATA cipr(1), cipi(1), cipr(2), cipi(2), cipr(3), cipi(3), cipr(4), &
    cipi(4)/1.0D0, 0.0D0, 0.0D0, 1.0D0, -1.0D0, 0.0D0, 0.0D0, &
    -1.0D0/
  DATA hpi, aic/1.57079632679489662D+00, 1.265512123484645396D+00/
  !***FIRST EXECUTABLE STATEMENT  ZUNI2
  Nz = 0
  nd = N
  Nlast = 0
  !-----------------------------------------------------------------------
  !     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
  !     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
  !     EXP(ALIM)=EXP(ELIM)*TOL
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
  !-----------------------------------------------------------------------
  !     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
  !-----------------------------------------------------------------------
  znr = Zi
  zni = -Zr
  zbr = Zr
  zbi = Zi
  cidi = -coner
  inu = Fnu
  ang = hpi*(Fnu-inu)
  c2r = COS(ang)
  c2i = SIN(ang)
  car = c2r
  sar = c2i
  in = inu + N - 1
  in = MOD(in,4) + 1
  str = c2r*cipr(in) - c2i*cipi(in)
  c2i = c2r*cipi(in) + c2i*cipr(in)
  c2r = str
  IF ( Zi<=0.0D0 ) THEN
    znr = -znr
    zbi = -zbi
    cidi = -cidi
    c2i = -c2i
  ENDIF
  !-----------------------------------------------------------------------
  !     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
  !-----------------------------------------------------------------------
  fn = MAX(Fnu,1.0D0)
  CALL ZUNHJ(znr,zni,fn,1,Tol,phir,phii,argr,argi,zeta1r,zeta1i,zeta2r,&
    zeta2i,asumr,asumi,bsumr,bsumi)
  IF ( Kode==1 ) THEN
    s1r = -zeta1r + zeta2r
    s1i = -zeta1i + zeta2i
  ELSE
    str = zbr + zeta2r
    sti = zbi + zeta2i
    rast = fn/ZABS(str,sti)
    str = str*rast*rast
    sti = -sti*rast*rast
    s1r = -zeta1r + str
    s1i = -zeta1i + sti
  ENDIF
  rs1 = s1r
  IF ( ABS(rs1)>Elim ) THEN
    IF ( rs1>0.0D0 ) GOTO 400
    Nz = N
    DO i = 1, N
      Yr(i) = zeror
      Yi(i) = zeroi
    ENDDO
    GOTO 99999
  ENDIF
  100  nn = MIN(2,nd)
  DO i = 1, nn
    fn = Fnu + (nd-i)
    CALL ZUNHJ(znr,zni,fn,0,Tol,phir,phii,argr,argi,zeta1r,zeta1i,zeta2r,&
      zeta2i,asumr,asumi,bsumr,bsumi)
    IF ( Kode==1 ) THEN
      s1r = -zeta1r + zeta2r
      s1i = -zeta1i + zeta2i
    ELSE
      str = zbr + zeta2r
      sti = zbi + zeta2i
      rast = fn/ZABS(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = -zeta1r + str
      s1i = -zeta1i + sti + ABS(Zi)
    ENDIF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = s1r
    IF ( ABS(rs1)>Elim ) GOTO 300
    IF ( i==1 ) iflag = 2
    IF ( ABS(rs1)>=Alim ) THEN
      !-----------------------------------------------------------------------
      !     REFINE  TEST AND SCALE
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      aphi = ZABS(phir,phii)
      aarg = ZABS(argr,argi)
      rs1 = rs1 + LOG(aphi) - 0.25D0*LOG(aarg) - aic
      IF ( ABS(rs1)>Elim ) GOTO 300
      IF ( i==1 ) iflag = 1
      IF ( rs1>=0.0D0 ) THEN
        IF ( i==1 ) iflag = 3
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------
    !     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
    !     EXPONENT EXTREMES
    !-----------------------------------------------------------------------
    CALL ZAIRY(argr,argi,0,2,air,aii,nai,idum)
    CALL ZAIRY(argr,argi,1,2,dair,daii,ndai,idum)
    str = dair*bsumr - daii*bsumi
    sti = dair*bsumi + daii*bsumr
    str = str + (air*asumr-aii*asumi)
    sti = sti + (air*asumi+aii*asumr)
    s2r = phir*str - phii*sti
    s2i = phir*sti + phii*str
    str = EXP(s1r)*cssr(iflag)
    s1r = str*COS(s1i)
    s1i = str*SIN(s1i)
    str = s2r*s1r - s2i*s1i
    s2i = s2r*s1i + s2i*s1r
    s2r = str
    IF ( iflag==1 ) THEN
      CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
      IF ( nw/=0 ) GOTO 300
    ENDIF
    IF ( Zi<=0.0D0 ) s2i = -s2i
    str = s2r*c2r - s2i*c2i
    s2i = s2r*c2i + s2i*c2r
    s2r = str
    cyr(i) = s2r
    cyi(i) = s2i
    j = nd - i + 1
    Yr(j) = s2r*csrr(iflag)
    Yi(j) = s2i*csrr(iflag)
    str = -c2i*cidi
    c2i = c2r*cidi
    c2r = str
  ENDDO
  IF ( nd>2 ) THEN
    raz = 1.0D0/ZABS(Zr,Zi)
    str = Zr*raz
    sti = -Zi*raz
    rzr = (str+str)*raz
    rzi = (sti+sti)*raz
    bry(2) = 1.0D0/bry(1)
    bry(3) = D1MACH(2)
    s1r = cyr(1)
    s1i = cyi(1)
    s2r = cyr(2)
    s2i = cyi(2)
    c1r = csrr(iflag)
    ascle = bry(iflag)
    k = nd - 2
    fn = k
    DO i = 3, nd
      c2r = s2r
      c2i = s2i
      s2r = s1r + (Fnu+fn)*(rzr*c2r-rzi*c2i)
      s2i = s1i + (Fnu+fn)*(rzr*c2i+rzi*c2r)
      s1r = c2r
      s1i = c2i
      c2r = s2r*c1r
      c2i = s2i*c1r
      Yr(k) = c2r
      Yi(k) = c2i
      k = k - 1
      fn = fn - 1.0D0
      IF ( iflag<3 ) THEN
        str = ABS(c2r)
        sti = ABS(c2i)
        c2m = MAX(str,sti)
        IF ( c2m>ascle ) THEN
          iflag = iflag + 1
          ascle = bry(iflag)
          s1r = s1r*c1r
          s1i = s1i*c1r
          s2r = c2r
          s2i = c2i
          s1r = s1r*cssr(iflag)
          s1i = s1i*cssr(iflag)
          s2r = s2r*cssr(iflag)
          s2i = s2i*cssr(iflag)
          c1r = csrr(iflag)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  200  RETURN
  300 CONTINUE
  IF ( rs1<=0.0D0 ) THEN
    !-----------------------------------------------------------------------
    !     SET UNDERFLOW AND UPDATE PARAMETERS
    !-----------------------------------------------------------------------
    Yr(nd) = zeror
    Yi(nd) = zeroi
    Nz = Nz + 1
    nd = nd - 1
    IF ( nd==0 ) GOTO 200
    CALL ZUOIK(Zr,Zi,Fnu,Kode,1,nd,Yr,Yi,nuf,Tol,Elim,Alim)
    IF ( nuf>=0 ) THEN
      nd = nd - nuf
      Nz = Nz + nuf
      IF ( nd==0 ) GOTO 200
      fn = Fnu + (nd-1)
      IF ( fn<Fnul ) THEN
        Nlast = nd
        RETURN
      ELSE
        !      FN = CIDI
        !      J = NUF + 1
        !      K = MOD(J,4) + 1
        !      S1R = CIPR(K)
        !      S1I = CIPI(K)
        !      IF (FN.LT.0.0D0) S1I = -S1I
        !      STR = C2R*S1R - C2I*S1I
        !      C2I = C2R*S1I + C2I*S1R
        !      C2R = STR
        in = inu + nd - 1
        in = MOD(in,4) + 1
        c2r = car*cipr(in) - sar*cipi(in)
        c2i = car*cipi(in) + sar*cipr(in)
        IF ( Zi<=0.0D0 ) c2i = -c2i
        GOTO 100
      ENDIF
    ENDIF
  ENDIF
  400  Nz = -1
  RETURN
  99999 CONTINUE
  END SUBROUTINE ZUNI2
