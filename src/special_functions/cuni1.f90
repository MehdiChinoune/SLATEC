!DECK CUNI1
SUBROUTINE CUNI1(Z,Fnu,Kode,N,Y,Nz,Nlast,Fnul,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CUNI1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CUNI1-A, ZUNI1-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
  !     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
  !
  !     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
  !     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
  !     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
  !     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
  !     Y(I)=CZERO FOR I=NLAST+1,N
  !
  !***SEE ALSO  CBESI, CBESK
  !***ROUTINES CALLED  CUCHK, CUNIK, CUOIK, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CUNI1
  COMPLEX cfn, cone, crsc, cscl, csr, css, cwrk, czero, c1, c2, &
    phi, rz, sum, s1, s2, Y, Z, zeta1, zeta2, cy
  REAL Alim, aphi, ascle, bry, c2i, c2m, c2r, Elim, fn, Fnu, &
    Fnul, rs1, Tol, yy, R1MACH
  INTEGER i, iflag, init, k, Kode, m, N, nd, Nlast, nn, nuf, nw, &
    Nz
  DIMENSION bry(3), Y(N), cwrk(16), css(3), csr(3), cy(2)
  DATA czero, cone/(0.0E0,0.0E0), (1.0E0,0.0E0)/
  !***FIRST EXECUTABLE STATEMENT  CUNI1
  Nz = 0
  nd = N
  Nlast = 0
  !-----------------------------------------------------------------------
  !     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
  !     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
  !     EXP(ALIM)=EXP(ELIM)*TOL
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
  !-----------------------------------------------------------------------
  !     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
  !-----------------------------------------------------------------------
  fn = MAX(Fnu,1.0E0)
  init = 0
  CALL CUNIK(Z,fn,1,1,Tol,init,phi,zeta1,zeta2,sum,cwrk)
  IF ( Kode==1 ) THEN
    s1 = -zeta1 + zeta2
  ELSE
    cfn = CMPLX(fn,0.0E0)
    s1 = -zeta1 + cfn*(cfn/(Z+zeta2))
  ENDIF
  rs1 = REAL(s1)
  IF ( ABS(rs1)>Elim ) THEN
    IF ( rs1>0.0E0 ) GOTO 400
    Nz = N
    DO i = 1, N
      Y(i) = czero
    ENDDO
    RETURN
  ENDIF
  100  nn = MIN(2,nd)
  DO i = 1, nn
    fn = Fnu + (nd-i)
    init = 0
    CALL CUNIK(Z,fn,1,0,Tol,init,phi,zeta1,zeta2,sum,cwrk)
    IF ( Kode==1 ) THEN
      s1 = -zeta1 + zeta2
    ELSE
      cfn = CMPLX(fn,0.0E0)
      yy = AIMAG(Z)
      s1 = -zeta1 + cfn*(cfn/(Z+zeta2)) + CMPLX(0.0E0,yy)
    ENDIF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1)
    IF ( ABS(rs1)>Elim ) GOTO 300
    IF ( i==1 ) iflag = 2
    IF ( ABS(rs1)>=Alim ) THEN
      !-----------------------------------------------------------------------
      !     REFINE  TEST AND SCALE
      !-----------------------------------------------------------------------
      aphi = ABS(phi)
      rs1 = rs1 + ALOG(aphi)
      IF ( ABS(rs1)>Elim ) GOTO 300
      IF ( i==1 ) iflag = 1
      IF ( rs1>=0.0E0 ) THEN
        IF ( i==1 ) iflag = 3
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------
    !     SCALE S1 IF ABS(S1).LT.ASCLE
    !-----------------------------------------------------------------------
    s2 = phi*sum
    c2r = REAL(s1)
    c2i = AIMAG(s1)
    c2m = EXP(c2r)*REAL(css(iflag))
    s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
    s2 = s2*s1
    IF ( iflag==1 ) THEN
      CALL CUCHK(s2,nw,bry(1),Tol)
      IF ( nw/=0 ) GOTO 300
    ENDIF
    m = nd - i + 1
    cy(i) = s2
    Y(m) = s2*csr(iflag)
  ENDDO
  IF ( nd>2 ) THEN
    rz = CMPLX(2.0E0,0.0E0)/Z
    bry(2) = 1.0E0/bry(1)
    bry(3) = R1MACH(2)
    s1 = cy(1)
    s2 = cy(2)
    c1 = csr(iflag)
    ascle = bry(iflag)
    k = nd - 2
    fn = k
    DO i = 3, nd
      c2 = s2
      s2 = s1 + CMPLX(Fnu+fn,0.0E0)*rz*s2
      s1 = c2
      c2 = s2*c1
      Y(k) = c2
      k = k - 1
      fn = fn - 1.0E0
      IF ( iflag<3 ) THEN
        c2r = REAL(c2)
        c2i = AIMAG(c2)
        c2r = ABS(c2r)
        c2i = ABS(c2i)
        c2m = MAX(c2r,c2i)
        IF ( c2m>ascle ) THEN
          iflag = iflag + 1
          ascle = bry(iflag)
          s1 = s1*c1
          s2 = c2
          s1 = s1*css(iflag)
          s2 = s2*css(iflag)
          c1 = csr(iflag)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  200  RETURN
  !-----------------------------------------------------------------------
  !     SET UNDERFLOW AND UPDATE PARAMETERS
  !-----------------------------------------------------------------------
  300 CONTINUE
  IF ( rs1<=0.0E0 ) THEN
    Y(nd) = czero
    Nz = Nz + 1
    nd = nd - 1
    IF ( nd==0 ) GOTO 200
    CALL CUOIK(Z,Fnu,Kode,1,nd,Y,nuf,Tol,Elim,Alim)
    IF ( nuf>=0 ) THEN
      nd = nd - nuf
      Nz = Nz + nuf
      IF ( nd==0 ) GOTO 200
      fn = Fnu + (nd-1)
      IF ( fn>=Fnul ) GOTO 100
      Nlast = nd
      RETURN
    ENDIF
  ENDIF
  400  Nz = -1
  RETURN
END SUBROUTINE CUNI1