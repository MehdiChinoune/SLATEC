!** ZKSCL
PURE SUBROUTINE ZKSCL(Zr,Fnu,N,Y,Nz,Rz,Ascle,Tol,Elim)
  !> Subsidiary to ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CKSCL-A, ZKSCL-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
  !     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
  !     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
  !
  !***
  ! **See also:**  ZBESK
  !***
  ! **Routines called:**  ZABS, ZLOG, ZUCHK

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)

  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Ascle, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Rz, Zr
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, ic, k, kk, nn, nw
  COMPLEX(DP) :: ck, cs, cy(2), s1, s2, zd, celm
  REAL(DP) :: aa, acs, as, csi, csr, fn, xx, zri, elm, alas, helim
  !* FIRST EXECUTABLE STATEMENT  ZKSCL
  Nz = 0
  ic = 0
  xx = REAL(Zr,DP)
  nn = MIN(2,N)
  DO i = 1, nn
    s1 = Y(i)
    cy(i) = s1
    as = ABS(s1)
    acs = -xx + LOG(as)
    Nz = Nz + 1
    Y(i) = (0._DP,0._DP)
    IF( acs>=(-Elim) ) THEN
      cs = -Zr + LOG(s1)
      csr = REAL(cs,DP)
      csi = AIMAG(cs)
      aa = EXP(csr)/Tol
      cs = CMPLX(aa,0._DP,DP)*CMPLX(COS(csi),SIN(csi),DP)
      CALL ZUCHK(cs,nw,Ascle,Tol)
      IF( nw==0 ) THEN
        Y(i) = cs
        Nz = Nz - 1
        ic = i
      END IF
    END IF
  END DO
  IF( N==1 ) RETURN
  IF( ic<=1 ) THEN
    Y(1) = (0._DP,0._DP)
    Nz = 2
  END IF
  IF( N==2 ) RETURN
  IF( Nz==0 ) RETURN
  fn = Fnu + 1._DP
  ck = CMPLX(fn,0._DP,DP)*Rz
  s1 = cy(1)
  s2 = cy(2)
  helim = 0.5_DP*Elim
  elm = EXP(-Elim)
  celm = CMPLX(elm,0._DP,DP)
  zri = AIMAG(Zr)
  zd = Zr
  !
  !     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
  !     S2 GETS LARGER THAN EXP(ELIM/2)
  !
  DO i = 3, N
    kk = i
    cs = s2
    s2 = ck*s2 + s1
    s1 = cs
    ck = ck + Rz
    as = ABS(s2)
    alas = LOG(as)
    acs = -xx + alas
    Nz = Nz + 1
    Y(i) = (0._DP,0._DP)
    IF( acs>=(-Elim) ) THEN
      cs = -zd + LOG(s2)
      csr = REAL(cs,DP)
      csi = AIMAG(cs)
      aa = EXP(csr)/Tol
      cs = CMPLX(aa,0._DP,DP)*CMPLX(COS(csi),SIN(csi),DP)
      CALL ZUCHK(cs,nw,Ascle,Tol)
      IF( nw==0 ) THEN
        Y(i) = cs
        Nz = Nz - 1
        IF( ic==(kk-1) ) GOTO 100
        ic = kk
        CYCLE
      END IF
    END IF
    IF( alas>=helim ) THEN
      xx = xx - Elim
      s1 = s1*celm
      s2 = s2*celm
      zd = CMPLX(xx,zri,DP)
    END IF
  END DO
  Nz = N
  IF( ic==N ) Nz = N - 1
  GOTO 200
  100  Nz = kk - 2
  200 CONTINUE
  DO k = 1, Nz
    Y(k) = (0._DP,0._DP)
  END DO
  !
END SUBROUTINE ZKSCL