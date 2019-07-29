!** ZUNI1
PURE SUBROUTINE ZUNI1(Z,Fnu,Kode,N,Y,Nz,Nlast,Fnul,Tol,Elim,Alim)
  !> Subsidiary to ZBESI and ZBESK
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
  !     EXPANSION FOR I(FNU,Z) IN -PI/3<=ARG Z<=PI/3.
  !
  !     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
  !     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
  !     NLAST/=0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
  !     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1<FNUL.
  !     Y(I)=CZERO FOR I=NLAST+1,N
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZUCHK, ZUNIK, ZUOIK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_dp, huge_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nlast, Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Fnul, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, iflag, init, k, m, nd, nn, nuf, nw
  COMPLEX(DP) :: cfn, crsc, cscl, csr(3), css(3), cwrk(16), c1, c2, &
    phi, rz, summ, s1, s2, zeta1, zeta2, cy(2)
  REAL(DP) :: aphi, ascle, bry(3), c2i, c2m, c2r, fn, rs1, yy
  !* FIRST EXECUTABLE STATEMENT  ZUNI1
  Nz = 0
  nd = N
  Nlast = 0
  !-----------------------------------------------------------------------
  !     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
  !     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
  !     EXP(ALIM)=EXP(ELIM)*TOL
  !-----------------------------------------------------------------------
  cscl = CMPLX(1._DP/Tol,0._DP,DP)
  crsc = CMPLX(Tol,0._DP,DP)
  css(1) = cscl
  css(2) = (1._DP,0._DP)
  css(3) = crsc
  csr(1) = crsc
  csr(2) = (1._DP,0._DP)
  csr(3) = cscl
  bry(1) = 1.E+3_DP*tiny_dp/Tol
  !-----------------------------------------------------------------------
  !     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
  !-----------------------------------------------------------------------
  fn = MAX(Fnu,1._DP)
  init = 0
  CALL ZUNIK(Z,fn,1,1,Tol,init,phi,zeta1,zeta2,summ,cwrk)
  IF( Kode==1 ) THEN
    s1 = -zeta1 + zeta2
  ELSE
    cfn = CMPLX(fn,0._DP,DP)
    s1 = -zeta1 + cfn*(cfn/(Z+zeta2))
  END IF
  rs1 = REAL(s1,DP)
  IF( ABS(rs1)>Elim ) THEN
    IF( rs1>0._DP ) GOTO 400
    Nz = N
    DO i = 1, N
      Y(i) = (0._DP,0._DP)
    END DO
    RETURN
  END IF
  100  nn = MIN(2,nd)
  DO i = 1, nn
    fn = Fnu + (nd-i)
    init = 0
    CALL ZUNIK(Z,fn,1,0,Tol,init,phi,zeta1,zeta2,summ,cwrk)
    IF( Kode==1 ) THEN
      s1 = -zeta1 + zeta2
    ELSE
      cfn = CMPLX(fn,0._DP,DP)
      yy = AIMAG(Z)
      s1 = -zeta1 + cfn*(cfn/(Z+zeta2)) + CMPLX(0._DP,yy,DP)
    END IF
    !-----------------------------------------------------------------------
    !     TEST FOR UNDERFLOW AND OVERFLOW
    !-----------------------------------------------------------------------
    rs1 = REAL(s1,DP)
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
      IF( rs1>=0._DP ) THEN
        IF( i==1 ) iflag = 3
      END IF
    END IF
    !-----------------------------------------------------------------------
    !     SCALE S1 IF ABS(S1)<ASCLE
    !-----------------------------------------------------------------------
    s2 = phi*summ
    c2r = REAL(s1,DP)
    c2i = AIMAG(s1)
    c2m = EXP(c2r)*REAL(css(iflag),DP)
    s1 = CMPLX(c2m,0._DP,DP)*CMPLX(COS(c2i),SIN(c2i),DP)
    s2 = s2*s1
    IF( iflag==1 ) THEN
      CALL ZUCHK(s2,nw,bry(1),Tol)
      IF( nw/=0 ) GOTO 300
    END IF
    m = nd - i + 1
    cy(i) = s2
    Y(m) = s2*csr(iflag)
  END DO
  IF( nd>2 ) THEN
    rz = CMPLX(2._DP,0._DP,DP)/Z
    bry(2) = 1._DP/bry(1)
    bry(3) = huge_dp
    s1 = cy(1)
    s2 = cy(2)
    c1 = csr(iflag)
    ascle = bry(iflag)
    k = nd - 2
    fn = k
    DO i = 3, nd
      c2 = s2
      s2 = s1 + CMPLX(Fnu+fn,0._DP,DP)*rz*s2
      s1 = c2
      c2 = s2*c1
      Y(k) = c2
      k = k - 1
      fn = fn - 1._DP
      IF( iflag<3 ) THEN
        c2r = REAL(c2,DP)
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
  IF( rs1<=0._DP ) THEN
    Y(nd) = (0._DP,0._DP)
    Nz = Nz + 1
    nd = nd - 1
    IF( nd==0 ) GOTO 200
    CALL ZUOIK(Z,Fnu,Kode,1,nd,Y,nuf,Tol,Elim,Alim)
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
  !
  RETURN
END SUBROUTINE ZUNI1