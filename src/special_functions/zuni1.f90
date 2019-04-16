!** ZUNI1
SUBROUTINE ZUNI1(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Nlast,Fnul,Tol,Elim,Alim)
  !>
  !***
  !  Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUNI1-A, ZUNI1-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
  !     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
  !
  !     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
  !     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
  !     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
  !     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
  !     Y(I)=CZERO FOR I=NLAST+1,N
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  !     COMPLEX CFN,CONE,CRSC,CSCL,CSR,CSS,CWRK,CZERO,C1,C2,PHI,RZ,SUM,S1,
  !    *S2,Y,Z,ZETA1,ZETA2
  INTEGER i, iflag, init, k, Kode, m, N, nd, Nlast, nn, nuf, nw, Nz
  REAL(8) :: Alim, aphi, ascle, bry(3), crsc, cscl, csrr(3), cssr(3), &
    cwrki(16), cwrkr(16), c1r, c2i, c2m, c2r, Elim, fn, Fnu, Fnul, phii, phir, &
    rast, rs1, rzi, rzr, sti, str, sumi, sumr, s1i, s1r, s2i, s2r, Tol, &
    Yi(N), Yr(N), zeta1i, zeta1r, zeta2i, zeta2r, Zi, Zr, cyr(2), cyi(2)
  REAL(8), PARAMETER :: zeror = 0.0D0, zeroi = 0.0D0, coner = 1.0D0
  !* FIRST EXECUTABLE STATEMENT  ZUNI1
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
  !     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
  !-----------------------------------------------------------------------
  fn = MAX(Fnu,1.0D0)
  init = 0
  CALL ZUNIK(Zr,Zi,fn,1,1,Tol,init,phir,phii,zeta1r,zeta1i,zeta2r,zeta2i,&
    sumr,sumi,cwrkr,cwrki)
  IF ( Kode==1 ) THEN
    s1r = -zeta1r + zeta2r
    s1i = -zeta1i + zeta2i
  ELSE
    str = Zr + zeta2r
    sti = Zi + zeta2i
    rast = fn/ZABS(str,sti)
    str = str*rast*rast
    sti = -sti*rast*rast
    s1r = -zeta1r + str
    s1i = -zeta1i + sti
  END IF
  rs1 = s1r
  IF ( ABS(rs1)>Elim ) THEN
    IF ( rs1>0.0D0 ) GOTO 400
    Nz = N
    DO i = 1, N
      Yr(i) = zeror
      Yi(i) = zeroi
    END DO
    RETURN
  END IF
  100  nn = MIN(2,nd)
  DO i = 1, nn
    fn = Fnu + (nd-i)
    init = 0
    CALL ZUNIK(Zr,Zi,fn,1,0,Tol,init,phir,phii,zeta1r,zeta1i,zeta2r,zeta2i,&
      sumr,sumi,cwrkr,cwrki)
    IF ( Kode==1 ) THEN
      s1r = -zeta1r + zeta2r
      s1i = -zeta1i + zeta2i
    ELSE
      str = Zr + zeta2r
      sti = Zi + zeta2i
      rast = fn/ZABS(str,sti)
      str = str*rast*rast
      sti = -sti*rast*rast
      s1r = -zeta1r + str
      s1i = -zeta1i + sti + Zi
    END IF
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
      aphi = ZABS(phir,phii)
      rs1 = rs1 + LOG(aphi)
      IF ( ABS(rs1)>Elim ) GOTO 300
      IF ( i==1 ) iflag = 1
      IF ( rs1>=0.0D0 ) THEN
        IF ( i==1 ) iflag = 3
      END IF
    END IF
    !-----------------------------------------------------------------------
    !     SCALE S1 IF ABS(S1).LT.ASCLE
    !-----------------------------------------------------------------------
    s2r = phir*sumr - phii*sumi
    s2i = phir*sumi + phii*sumr
    str = EXP(s1r)*cssr(iflag)
    s1r = str*COS(s1i)
    s1i = str*SIN(s1i)
    str = s2r*s1r - s2i*s1i
    s2i = s2r*s1i + s2i*s1r
    s2r = str
    IF ( iflag==1 ) THEN
      CALL ZUCHK(s2r,s2i,nw,bry(1),Tol)
      IF ( nw/=0 ) GOTO 300
    END IF
    cyr(i) = s2r
    cyi(i) = s2i
    m = nd - i + 1
    Yr(m) = s2r*csrr(iflag)
    Yi(m) = s2i*csrr(iflag)
  END DO
  IF ( nd>2 ) THEN
    rast = 1.0D0/ZABS(Zr,Zi)
    str = Zr*rast
    sti = -Zi*rast
    rzr = (str+str)*rast
    rzi = (sti+sti)*rast
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
        END IF
      END IF
    END DO
  END IF
  200  RETURN
  !-----------------------------------------------------------------------
  !     SET UNDERFLOW AND UPDATE PARAMETERS
  !-----------------------------------------------------------------------
  300 CONTINUE
  IF ( rs1<=0.0D0 ) THEN
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
      IF ( fn>=Fnul ) GOTO 100
      Nlast = nd
      RETURN
    END IF
  END IF
  400  Nz = -1
  RETURN
END SUBROUTINE ZUNI1
