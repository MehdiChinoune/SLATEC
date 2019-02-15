!DECK CASYI
SUBROUTINE CASYI(Z,Fnu,Kode,N,Y,Nz,Rl,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CASYI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CASYI-A, ZASYI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
  !     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
  !     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
  !
  !***SEE ALSO  CBESI, CBESK
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CASYI
  COMPLEX ak1, ck, cone, cs1, cs2, cz, czero, dk, ez, p1, rz, &
    s2, Y, Z
  REAL aa, acz, aez, ak, Alim, arg, arm, atol, az, bb, bk, dfnu, &
    dnu2, Elim, fdn, Fnu, pi, Rl, rtpi, rtr1, s, sgn, sqk, &
    Tol, x, yy, R1MACH
  INTEGER i, ib, il, inu, j, jl, k, Kode, koded, m, N, nn, Nz
  DIMENSION Y(N)
  DATA pi, rtpi/3.14159265358979324E0, 0.159154943091895336E0/
  DATA czero, cone/(0.0E0,0.0E0), (1.0E0,0.0E0)/
  !***FIRST EXECUTABLE STATEMENT  CASYI
  Nz = 0
  az = ABS(Z)
  x = REAL(Z)
  arm = 1.0E+3*R1MACH(1)
  rtr1 = SQRT(arm)
  il = MIN(2,N)
  dfnu = Fnu + (N-il)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  ak1 = CMPLX(rtpi,0.0E0)/Z
  ak1 = CSQRT(ak1)
  cz = Z
  IF ( Kode==2 ) cz = Z - CMPLX(x,0.0E0)
  acz = REAL(cz)
  IF ( ABS(acz)>Elim ) THEN
    Nz = -1
    RETURN
  ELSE
    dnu2 = dfnu + dfnu
    koded = 1
    IF ( (ABS(acz)<=Alim).OR.(N<=2) ) THEN
      koded = 0
      ak1 = ak1*CEXP(cz)
    ENDIF
    fdn = 0.0E0
    IF ( dnu2>rtr1 ) fdn = dnu2*dnu2
    ez = Z*CMPLX(8.0E0,0.0E0)
    !-----------------------------------------------------------------------
    !     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
    !     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
    !     EXPANSION FOR THE IMAGINARY PART.
    !-----------------------------------------------------------------------
    aez = 8.0E0*az
    s = Tol/aez
    jl = Rl + Rl + 2
    yy = AIMAG(Z)
    p1 = czero
    IF ( yy/=0.0E0 ) THEN
      !-----------------------------------------------------------------------
      !     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
      !     SIGNIFICANCE WHEN FNU OR N IS LARGE
      !-----------------------------------------------------------------------
      inu = Fnu
      arg = (Fnu-inu)*pi
      inu = inu + N - il
      ak = -SIN(arg)
      bk = COS(arg)
      IF ( yy<0.0E0 ) bk = -bk
      p1 = CMPLX(ak,bk)
      IF ( MOD(inu,2)==1 ) p1 = -p1
    ENDIF
    DO k = 1, il
      sqk = fdn - 1.0E0
      atol = s*ABS(sqk)
      sgn = 1.0E0
      cs1 = cone
      cs2 = cone
      ck = cone
      ak = 0.0E0
      aa = 1.0E0
      bb = aez
      dk = ez
      DO j = 1, jl
        ck = ck*CMPLX(sqk,0.0E0)/dk
        cs2 = cs2 + ck
        sgn = -sgn
        cs1 = cs1 + ck*CMPLX(sgn,0.0E0)
        dk = dk + ez
        aa = aa*ABS(sqk)/bb
        bb = bb + aez
        ak = ak + 8.0E0
        sqk = sqk - ak
        IF ( aa<=atol ) GOTO 20
      ENDDO
      GOTO 100
      20       s2 = cs1
      IF ( x+x<Elim ) s2 = s2 + p1*cs2*CEXP(-Z-Z)
      fdn = fdn + 8.0E0*dfnu + 4.0E0
      p1 = -p1
      m = N - il + k
      Y(m) = s2*ak1
    ENDDO
    IF ( N<=2 ) RETURN
    nn = N
    k = nn - 2
    ak = k
    rz = (cone+cone)/Z
    ib = 3
    DO i = ib, nn
      Y(k) = CMPLX(ak+Fnu,0.0E0)*rz*Y(k+1) + Y(k+2)
      ak = ak - 1.0E0
      k = k - 1
    ENDDO
    IF ( koded==0 ) RETURN
    ck = CEXP(cz)
    DO i = 1, nn
      Y(i) = Y(i)*ck
    ENDDO
    RETURN
  ENDIF
  100  Nz = -2
END SUBROUTINE CASYI
