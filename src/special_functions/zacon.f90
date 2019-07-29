!** ZACON
PURE SUBROUTINE ZACON(Z,Fnu,Kode,Mr,N,Y,Nz,Rl,Fnul,Tol,Elim,Alim)
  !> Subsidiary to ZBESH and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CACON-A, ZACON-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZACON APPLIES THE ANALYTIC CONTINUATION FORMULA
  !
  !         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
  !                 MP=PI*MR*CMPLX(0.0,1.0)
  !
  !     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
  !     HALF Z PLANE
  !
  !***
  ! **See also:**  ZBESH, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZBINU, ZBKNU, ZMLT, ZS1S2

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_dp, huge_dp
  USE IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE
  !
  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Fnul, Rl, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, inu, iuf, kflag, nn, nw
  REAL(DP) :: arg, ascle, as2, bscle, bry(3), cpn, c1i, c1m, c1r, fmr, sgn, spn, yy
  COMPLEX(DP) :: ck, cs, cscl, cscr, csgn, cspn, css(3), csr(3), c1, c2, &
    rz, sc1, sc2, st, s1, s2, zn, cy(2)
  REAL(DP), PARAMETER :: pi = 3.14159265358979324_DP
  REAL(DP), PARAMETER :: sqrt_huge = SQRT( huge_dp )
  !* FIRST EXECUTABLE STATEMENT  ZACON
  Nz = 0
  zn = -Z
  nn = N
  CALL ZBINU(zn,Fnu,Kode,nn,Y,nw,Rl,Fnul,Tol,Elim,Alim)
  IF( nw>=0 ) THEN
    !-----------------------------------------------------------------------
    !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
    !-----------------------------------------------------------------------
    nn = MIN(2,N)
    CALL ZBKNU(zn,Fnu,Kode,nn,cy,nw,Tol,Elim,Alim)
    IF( nw==0 ) THEN
      s1 = cy(1)
      fmr = Mr
      sgn = -SIGN(pi,fmr)
      csgn = CMPLX(0._DP,sgn,DP)
      IF( Kode/=1 ) THEN
        yy = -AIMAG(zn)
        cpn = COS(yy)
        spn = SIN(yy)
        csgn = csgn*CMPLX(cpn,spn,DP)
      END IF
      !-----------------------------------------------------------------------
      !     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
      !     WHEN FNU IS LARGE
      !-----------------------------------------------------------------------
      inu = INT( Fnu )
      arg = (Fnu-inu)*sgn
      cpn = COS(arg)
      spn = SIN(arg)
      cspn = CMPLX(cpn,spn,DP)
      IF( MOD(inu,2)==1 ) cspn = -cspn
      iuf = 0
      c1 = s1
      c2 = Y(1)
      ascle = 1.E+3_DP*tiny_dp/Tol
      IF( Kode/=1 ) THEN
        CALL ZS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
        Nz = Nz + nw
        sc1 = c1
      END IF
      Y(1) = cspn*c1 + csgn*c2
      IF( N==1 ) RETURN
      cspn = -cspn
      s2 = cy(2)
      c1 = s2
      c2 = Y(2)
      IF( Kode/=1 ) THEN
        CALL ZS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
        Nz = Nz + nw
        sc2 = c1
      END IF
      Y(2) = cspn*c1 + csgn*c2
      IF( N==2 ) RETURN
      cspn = -cspn
      rz = CMPLX(2._DP,0._DP,DP)/zn
      ck = CMPLX(Fnu+1._DP,0._DP,DP)*rz
      !-----------------------------------------------------------------------
      !     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
      !-----------------------------------------------------------------------
      cscl = CMPLX(1._DP/Tol,0._DP,DP)
      cscr = CMPLX(Tol,0._DP,DP)
      css(1) = cscl
      css(2) = (1._DP,0._DP)
      css(3) = cscr
      csr(1) = cscr
      csr(2) = (1._DP,0._DP)
      csr(3) = cscl
      bry(1) = ascle
      bry(2) = 1._DP/ascle
      bry(3) = huge_dp
      as2 = ABS(s2)
      IF( .NOT. IEEE_IS_FINITE(as2) ) as2 = ABS(s2/sqrt_huge) * sqrt_huge
      kflag = 2
      IF( as2<=bry(1) ) THEN
        kflag = 1
      ELSEIF( as2>=bry(2) ) THEN
        kflag = 3
      END IF
      bscle = bry(kflag)
      s1 = s1*css(kflag)
      s2 = s2*css(kflag)
      cs = csr(kflag)
      DO i = 3, N
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        c1 = s2*cs
        st = c1
        c2 = Y(i)
        IF( Kode/=1 ) THEN
          IF( iuf>=0 ) THEN
            CALL ZS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
            Nz = Nz + nw
            sc1 = sc2
            sc2 = c1
            IF( iuf==3 ) THEN
              iuf = -4
              s1 = sc1*css(kflag)
              s2 = sc2*css(kflag)
              st = sc2
            END IF
          END IF
        END IF
        Y(i) = cspn*c1 + csgn*c2
        ck = ck + rz
        cspn = -cspn
        IF( kflag<3 ) THEN
          c1r = REAL(c1,DP)
          c1i = AIMAG(c1)
          c1r = ABS(c1r)
          c1i = ABS(c1i)
          c1m = MAX(c1r,c1i)
          IF( c1m>bscle ) THEN
            kflag = kflag + 1
            bscle = bry(kflag)
            s1 = s1*cs
            s2 = st
            s1 = s1*css(kflag)
            s2 = s2*css(kflag)
            cs = csr(kflag)
          END IF
        END IF
      END DO
      RETURN
    END IF
  END IF
  Nz = -1
  IF( nw==(-2) ) Nz = -2
  !
END SUBROUTINE ZACON