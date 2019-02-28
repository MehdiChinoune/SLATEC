!DECK ZSERI
SUBROUTINE ZSERI(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ZSERI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CSERI-A, ZSERI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
  !     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
  !     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
  !     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
  !     CONDITION ABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
  !     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
  !
  !***SEE ALSO  ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, DGAMLN, ZABS, ZDIV, ZLOG, ZMLT, ZUCHK
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZLOG to EXTERNAL statement.  (RWC)
  !***END PROLOGUE  ZSERI
  !     COMPLEX AK1,CK,COEF,CONE,CRSC,CSCL,CZ,CZERO,HZ,RZ,S1,S2,Y,Z
  REAL(8) :: aa, acz, ak, ak1i, ak1r, Alim, arm, ascle, atol, &
    az, cki, ckr, coefi, coefr, conei, coner, crscr, &
    czi, czr, dfnu, Elim, Fnu, fnup, hzi, hzr, raz, &
    rs, rtr1, rzi, rzr, s, ss, sti, str, s1i, s1r, &
    s2i, s2r, Tol, Yi, Yr, wi, wr, zeroi, zeror, &
    Zi, Zr, DGAMLN, D1MACH, ZABS
  INTEGER i, ib, idum, iflag, il, k, Kode, l, m, N, nn, Nz, nw
  DIMENSION Yr(N), Yi(N), wr(2), wi(2)
  EXTERNAL ZABS, ZLOG
  DATA zeror, zeroi, coner, conei/0.0D0, 0.0D0, 1.0D0, 0.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZSERI
  Nz = 0
  az = ZABS(Zr,Zi)
  IF ( az==0.0D0 ) GOTO 500
  arm = 1.0D+3*D1MACH(1)
  rtr1 = SQRT(arm)
  crscr = 1.0D0
  iflag = 0
  IF ( az<arm ) THEN
    Nz = N
    IF ( Fnu==0.0D0 ) Nz = Nz - 1
    GOTO 500
  ELSE
    hzr = 0.5D0*Zr
    hzi = 0.5D0*Zi
    czr = zeror
    czi = zeroi
    IF ( az>rtr1 ) CALL ZMLT(hzr,hzi,hzr,hzi,czr,czi)
    acz = ZABS(czr,czi)
    nn = N
    CALL ZLOG(hzr,hzi,ckr,cki,idum)
  ENDIF
  100  dfnu = Fnu + (nn-1)
  fnup = dfnu + 1.0D0
  !-----------------------------------------------------------------------
  !     UNDERFLOW TEST
  !-----------------------------------------------------------------------
  ak1r = ckr*dfnu
  ak1i = cki*dfnu
  ak = DGAMLN(fnup,idum)
  ak1r = ak1r - ak
  IF ( Kode==2 ) ak1r = ak1r - Zr
  IF ( ak1r>(-Elim) ) THEN
    IF ( ak1r<=(-Alim) ) THEN
      iflag = 1
      ss = 1.0D0/Tol
      crscr = Tol
      ascle = arm*ss
    ENDIF
    aa = EXP(ak1r)
    IF ( iflag==1 ) aa = aa*ss
    coefr = aa*COS(ak1i)
    coefi = aa*SIN(ak1i)
    atol = Tol*acz/fnup
    il = MIN(2,nn)
    DO i = 1, il
      dfnu = Fnu + (nn-i)
      fnup = dfnu + 1.0D0
      s1r = coner
      s1i = conei
      IF ( acz>=Tol*fnup ) THEN
        ak1r = coner
        ak1i = conei
        ak = fnup + 2.0D0
        s = fnup
        aa = 2.0D0
        DO
          rs = 1.0D0/s
          str = ak1r*czr - ak1i*czi
          sti = ak1r*czi + ak1i*czr
          ak1r = str*rs
          ak1i = sti*rs
          s1r = s1r + ak1r
          s1i = s1i + ak1i
          s = s + ak
          ak = ak + 2.0D0
          aa = aa*acz*rs
          IF ( aa<=atol ) EXIT
        ENDDO
      ENDIF
      s2r = s1r*coefr - s1i*coefi
      s2i = s1r*coefi + s1i*coefr
      wr(i) = s2r
      wi(i) = s2i
      IF ( iflag/=0 ) THEN
        CALL ZUCHK(s2r,s2i,nw,ascle,Tol)
        IF ( nw/=0 ) GOTO 200
      ENDIF
      m = nn - i + 1
      Yr(m) = s2r*crscr
      Yi(m) = s2i*crscr
      IF ( i/=il ) THEN
        CALL ZDIV(coefr,coefi,hzr,hzi,str,sti)
        coefr = str*dfnu
        coefi = sti*dfnu
      ENDIF
    ENDDO
    IF ( nn<=2 ) RETURN
    k = nn - 2
    ak = k
    raz = 1.0D0/az
    str = Zr*raz
    sti = -Zi*raz
    rzr = (str+str)*raz
    rzi = (sti+sti)*raz
    IF ( iflag==1 ) THEN
      !-----------------------------------------------------------------------
      !     RECUR BACKWARD WITH SCALED VALUES
      !-----------------------------------------------------------------------
      !-----------------------------------------------------------------------
      !     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
      !     UNDERFLOW LIMIT = ASCLE = D1MACH(1)*SS*1.0D+3
      !-----------------------------------------------------------------------
      s1r = wr(1)
      s1i = wi(1)
      s2r = wr(2)
      s2i = wi(2)
      DO l = 3, nn
        ckr = s2r
        cki = s2i
        s2r = s1r + (ak+Fnu)*(rzr*ckr-rzi*cki)
        s2i = s1i + (ak+Fnu)*(rzr*cki+rzi*ckr)
        s1r = ckr
        s1i = cki
        ckr = s2r*crscr
        cki = s2i*crscr
        Yr(k) = ckr
        Yi(k) = cki
        ak = ak - 1.0D0
        k = k - 1
        IF ( ZABS(ckr,cki)>ascle ) GOTO 400
      ENDDO
      RETURN
    ELSE
      ib = 3
      GOTO 300
    ENDIF
  ENDIF
  200  Nz = Nz + 1
  Yr(nn) = zeror
  Yi(nn) = zeroi
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
    Yr(k) = (ak+Fnu)*(rzr*Yr(k+1)-rzi*Yi(k+1)) + Yr(k+2)
    Yi(k) = (ak+Fnu)*(rzr*Yi(k+1)+rzi*Yr(k+1)) + Yi(k+2)
    ak = ak - 1.0D0
    k = k - 1
  ENDDO
  RETURN
  400  ib = l + 1
  IF ( ib>nn ) RETURN
  GOTO 300
  500  Yr(1) = zeror
  Yi(1) = zeroi
  IF ( Fnu==0.0D0 ) THEN
    Yr(1) = coner
    Yi(1) = conei
  ENDIF
  IF ( N==1 ) RETURN
  DO i = 2, N
    Yr(i) = zeror
    Yi(i) = zeroi
  ENDDO
  RETURN
END SUBROUTINE ZSERI
