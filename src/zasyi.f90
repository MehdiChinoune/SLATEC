!*==ZASYI.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZASYI
SUBROUTINE ZASYI(Zr,Zi,Fnu,Kode,N,Yr,Yi,Nz,Rl,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--ZASYI5
  !***BEGIN PROLOGUE  ZASYI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZBESI and ZBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CASYI-A, ZASYI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
  !     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
  !     REGION ABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
  !     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
  !
  !***SEE ALSO  ZBESI, ZBESK
  !***ROUTINES CALLED  D1MACH, ZABS, ZDIV, ZEXP, ZMLT, ZSQRT
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !   930122  Added ZEXP and ZSQRT to EXTERNAL statement.  (RWC)
  !***END PROLOGUE  ZASYI
  !     COMPLEX AK1,CK,CONE,CS1,CS2,CZ,CZERO,DK,EZ,P1,RZ,S2,Y,Z
  DOUBLE PRECISION aa , aez , ak , ak1i , ak1r , Alim , arg , arm , atol , &
    az , bb , bk , cki , ckr , conei , coner , cs1i , cs1r , &
    cs2i , cs2r , czi , czr , dfnu , dki , dkr , dnu2 , &
    Elim , ezi , ezr , fdn , Fnu , pi , p1i , p1r , raz , &
    Rl , rtpi , rtr1 , rzi , rzr , s , sgn , sqk , sti , &
    str , s2i , s2r , Tol , tzi , tzr , Yi , Yr , zeroi , &
    zeror , Zi , Zr , D1MACH , ZABS
  INTEGER i , ib , il , inu , j , jl , k , Kode , koded , m , N , nn , Nz
  DIMENSION Yr(N) , Yi(N)
  EXTERNAL ZABS , ZEXP , ZSQRT
  DATA pi , rtpi/3.14159265358979324D0 , 0.159154943091895336D0/
  DATA zeror , zeroi , coner , conei/0.0D0 , 0.0D0 , 1.0D0 , 0.0D0/
  !***FIRST EXECUTABLE STATEMENT  ZASYI
  Nz = 0
  az = ZABS(Zr,Zi)
  arm = 1.0D+3*D1MACH(1)
  rtr1 = SQRT(arm)
  il = MIN(2,N)
  dfnu = Fnu + (N-il)
  !-----------------------------------------------------------------------
  !     OVERFLOW TEST
  !-----------------------------------------------------------------------
  raz = 1.0D0/az
  str = Zr*raz
  sti = -Zi*raz
  ak1r = rtpi*str*raz
  ak1i = rtpi*sti*raz
  CALL ZSQRT(ak1r,ak1i,ak1r,ak1i)
  czr = Zr
  czi = Zi
  IF ( Kode==2 ) THEN
    czr = zeror
    czi = Zi
  ENDIF
  IF ( ABS(czr)>Elim ) THEN
    Nz = -1
    RETURN
  ELSE
    dnu2 = dfnu + dfnu
    koded = 1
    IF ( (ABS(czr)<=Alim).OR.(N<=2) ) THEN
      koded = 0
      CALL ZEXP(czr,czi,str,sti)
      CALL ZMLT(ak1r,ak1i,str,sti,ak1r,ak1i)
    ENDIF
    fdn = 0.0D0
    IF ( dnu2>rtr1 ) fdn = dnu2*dnu2
    ezr = Zr*8.0D0
    ezi = Zi*8.0D0
    !-----------------------------------------------------------------------
    !     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
    !     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
    !     EXPANSION FOR THE IMAGINARY PART.
    !-----------------------------------------------------------------------
    aez = 8.0D0*az
    s = Tol/aez
    jl = Rl + Rl + 2
    p1r = zeror
    p1i = zeroi
    IF ( Zi/=0.0D0 ) THEN
      !-----------------------------------------------------------------------
      !     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
      !     SIGNIFICANCE WHEN FNU OR N IS LARGE
      !-----------------------------------------------------------------------
      inu = Fnu
      arg = (Fnu-inu)*pi
      inu = inu + N - il
      ak = -SIN(arg)
      bk = COS(arg)
      IF ( Zi<0.0D0 ) bk = -bk
      p1r = ak
      p1i = bk
      IF ( MOD(inu,2)/=0 ) THEN
        p1r = -p1r
        p1i = -p1i
      ENDIF
    ENDIF
    DO k = 1 , il
      sqk = fdn - 1.0D0
      atol = s*ABS(sqk)
      sgn = 1.0D0
      cs1r = coner
      cs1i = conei
      cs2r = coner
      cs2i = conei
      ckr = coner
      cki = conei
      ak = 0.0D0
      aa = 1.0D0
      bb = aez
      dkr = ezr
      dki = ezi
      DO j = 1 , jl
        CALL ZDIV(ckr,cki,dkr,dki,str,sti)
        ckr = str*sqk
        cki = sti*sqk
        cs2r = cs2r + ckr
        cs2i = cs2i + cki
        sgn = -sgn
        cs1r = cs1r + ckr*sgn
        cs1i = cs1i + cki*sgn
        dkr = dkr + ezr
        dki = dki + ezi
        aa = aa*ABS(sqk)/bb
        bb = bb + aez
        ak = ak + 8.0D0
        sqk = sqk - ak
        IF ( aa<=atol ) GOTO 20
      ENDDO
      GOTO 100
      20       s2r = cs1r
      s2i = cs1i
      IF ( Zr+Zr<Elim ) THEN
        tzr = Zr + Zr
        tzi = Zi + Zi
        CALL ZEXP(-tzr,-tzi,str,sti)
        CALL ZMLT(str,sti,p1r,p1i,str,sti)
        CALL ZMLT(str,sti,cs2r,cs2i,str,sti)
        s2r = s2r + str
        s2i = s2i + sti
      ENDIF
      fdn = fdn + 8.0D0*dfnu + 4.0D0
      p1r = -p1r
      p1i = -p1i
      m = N - il + k
      Yr(m) = s2r*ak1r - s2i*ak1i
      Yi(m) = s2r*ak1i + s2i*ak1r
    ENDDO
    IF ( N<=2 ) RETURN
    nn = N
    k = nn - 2
    ak = k
    str = Zr*raz
    sti = -Zi*raz
    rzr = (str+str)*raz
    rzi = (sti+sti)*raz
    ib = 3
    DO i = ib , nn
      Yr(k) = (ak+Fnu)*(rzr*Yr(k+1)-rzi*Yi(k+1)) + Yr(k+2)
      Yi(k) = (ak+Fnu)*(rzr*Yi(k+1)+rzi*Yr(k+1)) + Yi(k+2)
      ak = ak - 1.0D0
      k = k - 1
    ENDDO
    IF ( koded==0 ) RETURN
    CALL ZEXP(czr,czi,ckr,cki)
    DO i = 1 , nn
      str = Yr(i)*ckr - Yi(i)*cki
      Yi(i) = Yr(i)*cki + Yi(i)*ckr
      Yr(i) = str
    ENDDO
    RETURN
  ENDIF
  100  Nz = -2
END SUBROUTINE ZASYI
