!** CSERI
SUBROUTINE CSERI(Z,Fnu,Kode,N,Y,Nz,Tol,Elim,Alim)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBESI and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CSERI-A, ZSERI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
  !     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
  !     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
  !     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
  !     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
  !     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
  !
  !***
  ! **See also:**  CBESI, CBESK
  !***
  ! **Routines called:**  CUCHK, GAMLN, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  
  INTEGER i, ib, idum, iflag, il, k, Kode, l, m, N, nn, nw, Nz
  COMPLEX ak1, ck, coef, crsc, cz, hz, rz, s1, s2, w(2), Y(N), Z
  REAL aa, acz, ak, Alim, arm, ascle, atol, az, dfnu, Elim, Fnu, &
    fnup, rak1, rs, rtr1, s, ss, Tol, x, GAMLN, R1MACH
  COMPLEX, PARAMETER :: czero = (0.0E0,0.0E0), cone = (1.0E0,0.0E0)
  !* FIRST EXECUTABLE STATEMENT  CSERI
  Nz = 0
  az = ABS(Z)
  IF ( az==0.0E0 ) GOTO 500
  x = REAL(Z)
  arm = 1.0E+3*R1MACH(1)
  rtr1 = SQRT(arm)
  crsc = CMPLX(1.0E0,0.0E0)
  iflag = 0
  IF ( az<arm ) THEN
    Nz = N
    IF ( Fnu==0.0E0 ) Nz = Nz - 1
    GOTO 500
  ELSE
    hz = Z*CMPLX(0.5E0,0.0E0)
    cz = czero
    IF ( az>rtr1 ) cz = hz*hz
    acz = ABS(cz)
    nn = N
    ck = CLOG(hz)
  ENDIF
  100  dfnu = Fnu + (nn-1)
  fnup = dfnu + 1.0E0
  !-----------------------------------------------------------------------
  !     UNDERFLOW TEST
  !-----------------------------------------------------------------------
  ak1 = ck*CMPLX(dfnu,0.0E0)
  ak = GAMLN(fnup,idum)
  ak1 = ak1 - CMPLX(ak,0.0E0)
  IF ( Kode==2 ) ak1 = ak1 - CMPLX(x,0.0E0)
  rak1 = REAL(ak1)
  IF ( rak1>(-Elim) ) THEN
    IF ( rak1<=(-Alim) ) THEN
      iflag = 1
      ss = 1.0E0/Tol
      crsc = CMPLX(Tol,0.0E0)
      ascle = arm*ss
    ENDIF
    ak = AIMAG(ak1)
    aa = EXP(rak1)
    IF ( iflag==1 ) aa = aa*ss
    coef = CMPLX(aa,0.0E0)*CMPLX(COS(ak),SIN(ak))
    atol = Tol*acz/fnup
    il = MIN(2,nn)
    DO i = 1, il
      dfnu = Fnu + (nn-i)
      fnup = dfnu + 1.0E0
      s1 = cone
      IF ( acz>=Tol*fnup ) THEN
        ak1 = cone
        ak = fnup + 2.0E0
        s = fnup
        aa = 2.0E0
        DO
          rs = 1.0E0/s
          ak1 = ak1*cz*CMPLX(rs,0.0E0)
          s1 = s1 + ak1
          s = s + ak
          ak = ak + 2.0E0
          aa = aa*acz*rs
          IF ( aa<=atol ) EXIT
        ENDDO
      ENDIF
      m = nn - i + 1
      s2 = s1*coef
      w(i) = s2
      IF ( iflag/=0 ) THEN
        CALL CUCHK(s2,nw,ascle,Tol)
        IF ( nw/=0 ) GOTO 200
      ENDIF
      Y(m) = s2*crsc
      IF ( i/=il ) coef = coef*CMPLX(dfnu,0.0E0)/hz
    ENDDO
    IF ( nn<=2 ) RETURN
    k = nn - 2
    ak = k
    rz = (cone+cone)/Z
    IF ( iflag==1 ) THEN
      !-----------------------------------------------------------------------
      !     RECUR BACKWARD WITH SCALED VALUES
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
      !     UNDERFLOW LIMIT = ASCLE = R1MACH(1)*CSCL*1.0E+3
      !-----------------------------------------------------------------------
      s1 = w(1)
      s2 = w(2)
      DO l = 3, nn
        ck = s2
        s2 = s1 + CMPLX(ak+Fnu,0.0E0)*rz*s2
        s1 = ck
        ck = s2*crsc
        Y(k) = ck
        ak = ak - 1.0E0
        k = k - 1
        IF ( ABS(ck)>ascle ) GOTO 400
      ENDDO
      RETURN
    ELSE
      ib = 3
      GOTO 300
    ENDIF
  ENDIF
  200  Nz = Nz + 1
  Y(nn) = czero
  IF ( acz>dfnu ) THEN
    !-----------------------------------------------------------------------
    !     RETURN WITH NZ.LT.0 IF ABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
    !     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
    !-----------------------------------------------------------------------
    Nz = -Nz
    RETURN
  ELSE
    nn = nn - 1
    IF ( nn==0 ) RETURN
    GOTO 100
  ENDIF
  300 CONTINUE
  DO i = ib, nn
    Y(k) = CMPLX(ak+Fnu,0.0E0)*rz*Y(k+1) + Y(k+2)
    ak = ak - 1.0E0
    k = k - 1
  ENDDO
  RETURN
  400  ib = l + 1
  IF ( ib>nn ) RETURN
  GOTO 300
  500  Y(1) = czero
  IF ( Fnu==0.0E0 ) Y(1) = cone
  IF ( N==1 ) RETURN
  DO i = 2, N
    Y(i) = czero
  ENDDO
  RETURN
END SUBROUTINE CSERI
