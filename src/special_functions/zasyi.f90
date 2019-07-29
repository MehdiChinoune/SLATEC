!** ZASYI
PURE SUBROUTINE ZASYI(Z,Fnu,Kode,N,Y,Nz,Rl,Tol,Elim,Alim)
  !> Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CASYI-A, ZASYI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z)>=0.0 BY
  !     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z)>MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
  !     NZ<0 INDICATES AN OVERFLOW ON KODE=1.
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
  USE service, ONLY : tiny_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Rl, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, ib, il, inu, j, jl, k, koded, m, nn
  COMPLEX(DP) :: ak1, ck, cs1, cs2, cz, dk, ez, p1, rz, s2
  REAL(DP) :: aa, acz, aez, ak, arg, arm, atol, az, bb, bk, dfnu, &
    dnu2, fdn, rtr1, s, sgn, sqk, x, yy
  REAL(DP), PARAMETER ::  pi = 3.14159265358979324_DP, rtpi = 0.159154943091895336_DP
  !* FIRST EXECUTABLE STATEMENT  ZASYI
  Nz = 0
  az = ABS(Z)
  x = REAL(Z,DP)
  arm = 1.E3_DP*tiny_dp
  rtr1 = SQRT(arm)
  il = MIN(2,N)
  dfnu = Fnu + (N-il)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  ak1 = CMPLX(rtpi,0._DP,DP)/Z
  ak1 = SQRT(ak1)
  cz = Z
  IF( Kode==2 ) cz = Z - CMPLX(x,0._DP,DP)
  acz = REAL(cz,DP)
  IF( ABS(acz)>Elim ) THEN
    Nz = -1
    RETURN
  ELSE
    dnu2 = dfnu + dfnu
    koded = 1
    IF( (ABS(acz)<=Alim) .OR. (N<=2) ) THEN
      koded = 0
      ak1 = ak1*EXP(cz)
    END IF
    fdn = 0._DP
    IF( dnu2>rtr1 ) fdn = dnu2*dnu2
    ez = Z*CMPLX(8._DP,0._DP,DP)
    !-----------------------------------------------------------------------
    !     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
    !     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
    !     EXPANSION FOR THE IMAGINARY PART.
    !-----------------------------------------------------------------------
    aez = 8._DP*az
    s = Tol/aez
    jl = INT( Rl + Rl ) + 2
    yy = AIMAG(Z)
    p1 = (0._DP,0._DP)
    IF( yy/=0._DP ) THEN
      !-----------------------------------------------------------------------
      !     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
      !     SIGNIFICANCE WHEN FNU OR N IS LARGE
      !-----------------------------------------------------------------------
      inu = INT( Fnu )
      arg = (Fnu-inu)*pi
      inu = inu + N - il
      ak = -SIN(arg)
      bk = COS(arg)
      IF( yy<0._DP ) bk = -bk
      p1 = CMPLX(ak,bk,DP)
      IF( MOD(inu,2)==1 ) p1 = -p1
    END IF
    DO k = 1, il
      sqk = fdn - 1._DP
      atol = s*ABS(sqk)
      sgn = 1._DP
      cs1 = (1._DP,0._DP)
      cs2 = (1._DP,0._DP)
      ck = (1._DP,0._DP)
      ak = 0._DP
      aa = 1._DP
      bb = aez
      dk = ez
      DO j = 1, jl
        ck = ck*CMPLX(sqk,0._DP,DP)/dk
        cs2 = cs2 + ck
        sgn = -sgn
        cs1 = cs1 + ck*CMPLX(sgn,0._DP,DP)
        dk = dk + ez
        aa = aa*ABS(sqk)/bb
        bb = bb + aez
        ak = ak + 8._DP
        sqk = sqk - ak
        IF( aa<=atol ) GOTO 20
      END DO
      GOTO 100
      20  s2 = cs1
      IF( x+x<Elim ) s2 = s2 + p1*cs2*EXP(-Z-Z)
      fdn = fdn + 8._DP*dfnu + 4._DP
      p1 = -p1
      m = N - il + k
      Y(m) = s2*ak1
    END DO
    IF( N<=2 ) RETURN
    nn = N
    k = nn - 2
    ak = k
    rz = (2._DP,0._DP)/Z
    ib = 3
    DO i = ib, nn
      Y(k) = CMPLX(ak+Fnu,0._DP,DP)*rz*Y(k+1) + Y(k+2)
      ak = ak - 1._DP
      k = k - 1
    END DO
    IF( koded==0 ) RETURN
    ck = EXP(cz)
    DO i = 1, nn
      Y(i) = Y(i)*ck
    END DO
    RETURN
  END IF
  100  Nz = -2
  !
END SUBROUTINE ZASYI