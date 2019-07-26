!** CASYI
PURE SUBROUTINE CASYI(Z,Fnu,Kode,N,Y,Nz,Rl,Tol,Elim,Alim)
  !> Subsidiary to CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CASYI-A, ZASYI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z)>=0.0 BY
  !     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z)>MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
  !     NZ<0 INDICATES AN OVERFLOW ON KODE=1.
  !
  !***
  ! **See also:**  CBESI, CBESK
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_sp
  !S
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alim, Elim, Fnu, Rl, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: i, ib, il, inu, j, jl, k, koded, m, nn
  COMPLEX(SP) :: ak1, ck, cs1, cs2, cz, dk, ez, p1, rz, s2
  REAL(SP) :: aa, acz, aez, ak, arg, arm, atol, az, bb, bk, dfnu, &
    dnu2, fdn, rtr1, s, sgn, sqk, x, yy
  REAL(SP), PARAMETER ::  pi = 3.14159265358979324_SP, rtpi = 0.159154943091895336_SP
  !* FIRST EXECUTABLE STATEMENT  CASYI
  Nz = 0
  az = ABS(Z)
  x = REAL(Z)
  arm = 1.E+3_SP*tiny_sp
  rtr1 = SQRT(arm)
  il = MIN(2,N)
  dfnu = Fnu + (N-il)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  ak1 = CMPLX(rtpi,0._SP,SP)/Z
  ak1 = SQRT(ak1)
  cz = Z
  IF( Kode==2 ) cz = Z - CMPLX(x,0._SP,SP)
  acz = REAL(cz)
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
    fdn = 0._SP
    IF( dnu2>rtr1 ) fdn = dnu2*dnu2
    ez = Z*CMPLX(8._SP,0._SP,SP)
    !-----------------------------------------------------------------------
    !     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
    !     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
    !     EXPANSION FOR THE IMAGINARY PART.
    !-----------------------------------------------------------------------
    aez = 8._SP*az
    s = Tol/aez
    jl = INT( Rl + Rl ) + 2
    yy = AIMAG(Z)
    p1 = (0._SP,0._SP)
    IF( yy/=0._SP ) THEN
      !-----------------------------------------------------------------------
      !     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
      !     SIGNIFICANCE WHEN FNU OR N IS LARGE
      !-----------------------------------------------------------------------
      inu = INT( Fnu )
      arg = (Fnu-inu)*pi
      inu = inu + N - il
      ak = -SIN(arg)
      bk = COS(arg)
      IF( yy<0._SP ) bk = -bk
      p1 = CMPLX(ak,bk,SP)
      IF( MOD(inu,2)==1 ) p1 = -p1
    END IF
    DO k = 1, il
      sqk = fdn - 1._SP
      atol = s*ABS(sqk)
      sgn = 1._SP
      cs1 = (1._SP,0._SP)
      cs2 = (1._SP,0._SP)
      ck = (1._SP,0._SP)
      ak = 0._SP
      aa = 1._SP
      bb = aez
      dk = ez
      DO j = 1, jl
        ck = ck*CMPLX(sqk,0._SP,SP)/dk
        cs2 = cs2 + ck
        sgn = -sgn
        cs1 = cs1 + ck*CMPLX(sgn,0._SP,SP)
        dk = dk + ez
        aa = aa*ABS(sqk)/bb
        bb = bb + aez
        ak = ak + 8._SP
        sqk = sqk - ak
        IF( aa<=atol ) GOTO 20
      END DO
      GOTO 100
      20  s2 = cs1
      IF( x+x<Elim ) s2 = s2 + p1*cs2*EXP(-Z-Z)
      fdn = fdn + 8._SP*dfnu + 4._SP
      p1 = -p1
      m = N - il + k
      Y(m) = s2*ak1
    END DO
    IF( N<=2 ) RETURN
    nn = N
    k = nn - 2
    ak = k
    rz = (2._SP,0._SP)/Z
    ib = 3
    DO i = ib, nn
      Y(k) = CMPLX(ak+Fnu,0._SP,SP)*rz*Y(k+1) + Y(k+2)
      ak = ak - 1._SP
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

END SUBROUTINE CASYI