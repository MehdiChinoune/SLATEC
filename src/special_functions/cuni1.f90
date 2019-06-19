!** CUNI1
SUBROUTINE CUNI1(Z,Fnu,Kode,N,Y,Nz,Nlast,Fnul,Tol,Elim,Alim)
  !> Subsidiary to CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CUNI1-A, ZUNI1-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
  !     EXPANSION FOR I(FNU,Z) IN -PI/3<=ARG Z<=PI/3.
  !
  !     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
  !     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
  !     NLAST/=0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
  !     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1<FNUL.
  !     Y(I)=CZERO FOR I=NLAST+1,N
  !
  !***
  ! **See also:**  CBESI, CBESK
  !***
  ! **Routines called:**  CUCHK, CUNIK, CUOIK, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  INTEGER :: i, iflag, init, k, Kode, m, N, nd, Nlast, nn, nuf, nw, Nz
  COMPLEX(SP) :: cfn, crsc, cscl, csr(3), css(3), cwrk(16), c1, c2, &
    phi, rz, summ, s1, s2, Y(N), Z, zeta1, zeta2, cy(2)
  REAL(SP) :: Alim, aphi, ascle, bry(3), c2i, c2m, c2r, Elim, fn, Fnu, &
    Fnul, rs1, Tol, yy
  COMPLEX(SP), PARAMETER :: czero = (0.0E0,0.0E0), cone = (1.0E0,0.0E0)
  !* FIRST EXECUTABLE STATEMENT  CUNI1
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
  CALL CUNIK(Z,fn,1,1,Tol,init,phi,zeta1,zeta2,summ,cwrk)
  IF( Kode==1 ) THEN
    s1 = -zeta1 + zeta2
  ELSE
    cfn = CMPLX(fn,0.0E0)
    s1 = -zeta1 + cfn*(cfn/(Z+zeta2))
  END IF
  rs1 = REAL(s1)
  IF( ABS(rs1)>Elim ) THEN
    IF( rs1>0.0E0 ) GOTO 400
    Nz = N
    DO i = 1, N
      Y(i) = czero
    END DO
    RETURN
  END IF
  100  nn = MIN(2,nd)
  DO i = 1, nn
    fn = Fnu + (nd-i)
    init = 0
    CALL CUNIK(Z,fn,1,0,Tol,init,phi,zeta1,zeta2,summ,cwrk)
    IF( Kode==1 ) THEN
      s1 = -zeta1 + zeta2
    ELSE
      cfn = CMPLX(fn,0.0E0)
      yy = AIMAG(Z)
      s1 = -zeta1 + cfn*(cfn/(Z+zeta2)) + CMPLX(0.0E0,yy)
    END IF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1)
    IF( ABS(rs1)>Elim ) GOTO 300
    IF( i==1 ) iflag = 2
    IF( ABS(rs1)>=Alim ) THEN
      !-----------------------------------------------------------------------
      !     REFINE  TEST AND SCALE
      !-----------------------------------------------------------------------
      aphi = ABS(phi)
      rs1 = rs1 + LOG(aphi)
      IF( ABS(rs1)>Elim ) GOTO 300
      IF( i==1 ) iflag = 1
      IF( rs1>=0.0E0 ) THEN
        IF( i==1 ) iflag = 3
      END IF
    END IF
    !-----------------------------------------------------------------------
    !     SCALE S1 IF ABS(S1)<ASCLE
    !-----------------------------------------------------------------------
    s2 = phi*summ
    c2r = REAL(s1)
    c2i = AIMAG(s1)
    c2m = EXP(c2r)*REAL(css(iflag))
    s1 = CMPLX(c2m,0.0E0)*CMPLX(COS(c2i),SIN(c2i))
    s2 = s2*s1
    IF( iflag==1 ) THEN
      CALL CUCHK(s2,nw,bry(1),Tol)
      IF( nw/=0 ) GOTO 300
    END IF
    m = nd - i + 1
    cy(i) = s2
    Y(m) = s2*csr(iflag)
  END DO
  IF( nd>2 ) THEN
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
      IF( iflag<3 ) THEN
        c2r = REAL(c2)
        c2i = AIMAG(c2)
        c2r = ABS(c2r)
        c2i = ABS(c2i)
        c2m = MAX(c2r,c2i)
        IF( c2m>ascle ) THEN
          iflag = iflag + 1
          ascle = bry(iflag)
          s1 = s1*c1
          s2 = c2
          s1 = s1*css(iflag)
          s2 = s2*css(iflag)
          c1 = csr(iflag)
        END IF
      END IF
    END DO
  END IF
  200  RETURN
  !-----------------------------------------------------------------------
  !     SET UNDERFLOW AND UPDATE PARAMETERS
  !-----------------------------------------------------------------------
  300 CONTINUE
  IF( rs1<=0.0E0 ) THEN
    Y(nd) = czero
    Nz = Nz + 1
    nd = nd - 1
    IF( nd==0 ) GOTO 200
    CALL CUOIK(Z,Fnu,Kode,1,nd,Y,nuf,Tol,Elim,Alim)
    IF( nuf>=0 ) THEN
      nd = nd - nuf
      Nz = Nz + nuf
      IF( nd==0 ) GOTO 200
      fn = Fnu + (nd-1)
      IF( fn>=Fnul ) GOTO 100
      Nlast = nd
      RETURN
    END IF
  END IF
  400  Nz = -1
  RETURN
END SUBROUTINE CUNI1
