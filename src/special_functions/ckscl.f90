!** CKSCL
SUBROUTINE CKSCL(Zr,Fnu,N,Y,Nz,Rz,Ascle,Tol,Elim)
  !>
  !  Subsidiary to CBKNU, CUNK1 and CUNK2
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
  ! **See also:**  CBKNU, CUNK1, CUNK2
  !***
  ! **Routines called:**  CUCHK

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER i, ic, k, kk, N, nn, nw, Nz
  COMPLEX(SP) ck, cs, cy(2), Rz, s1, s2, Y(N), Zr, zd, celm
  REAL(SP) aa, Ascle, acs, as, csi, csr, Elim, fn, Fnu, Tol, xx, zri, elm, alas, helim
  COMPLEX(SP), PARAMETER :: czero = (0.0E0,0.0E0)
  !* FIRST EXECUTABLE STATEMENT  CUCHK
  Nz = 0
  ic = 0
  xx = REAL(Zr)
  nn = MIN(2,N)
  DO i = 1, nn
    s1 = Y(i)
    cy(i) = s1
    as = ABS(s1)
    acs = -xx + LOG(as)
    Nz = Nz + 1
    Y(i) = czero
    IF ( acs>=(-Elim) ) THEN
      cs = -Zr + LOG(s1)
      csr = REAL(cs)
      csi = AIMAG(cs)
      aa = EXP(csr)/Tol
      cs = CMPLX(aa,0.0E0)*CMPLX(COS(csi),SIN(csi))
      CALL CUCHK(cs,nw,Ascle,Tol)
      IF ( nw==0 ) THEN
        Y(i) = cs
        Nz = Nz - 1
        ic = i
      END IF
    END IF
  END DO
  IF ( N==1 ) RETURN
  IF ( ic<=1 ) THEN
    Y(1) = czero
    Nz = 2
  END IF
  IF ( N==2 ) RETURN
  IF ( Nz==0 ) RETURN
  fn = Fnu + 1.0E0
  ck = CMPLX(fn,0.0E0)*Rz
  s1 = cy(1)
  s2 = cy(2)
  helim = 0.5E0*Elim
  elm = EXP(-Elim)
  celm = CMPLX(elm,0.0E0)
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
    Y(i) = czero
    IF ( acs>=(-Elim) ) THEN
      cs = -zd + LOG(s2)
      csr = REAL(cs)
      csi = AIMAG(cs)
      aa = EXP(csr)/Tol
      cs = CMPLX(aa,0.0E0)*CMPLX(COS(csi),SIN(csi))
      CALL CUCHK(cs,nw,Ascle,Tol)
      IF ( nw==0 ) THEN
        Y(i) = cs
        Nz = Nz - 1
        IF ( ic==(kk-1) ) GOTO 100
        ic = kk
        CYCLE
      END IF
    END IF
    IF ( alas>=helim ) THEN
      xx = xx - Elim
      s1 = s1*celm
      s2 = s2*celm
      zd = CMPLX(xx,zri)
    END IF
  END DO
  Nz = N
  IF ( ic==N ) Nz = N - 1
  GOTO 200
  100  Nz = kk - 2
  200 CONTINUE
  DO k = 1, Nz
    Y(k) = czero
  END DO
END SUBROUTINE CKSCL
