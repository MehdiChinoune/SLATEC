!*==ZACAI.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK ZACAI
SUBROUTINE ZACAI(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Rl,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--ZACAI5
  !***BEGIN PROLOGUE  ZACAI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to ZAIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CACAI-A, ZACAI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     ZACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
  !
  !         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
  !                 MP=PI*MR*CMPLX(0.0,1.0)
  !
  !     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
  !     HALF Z PLANE FOR USE WITH ZAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
  !     ZACAI IS THE SAME AS ZACON WITH THE PARTS FOR LARGER ORDERS AND
  !     RECURRENCE REMOVED. A RECURSIVE CALL TO ZACON CAN RESULT IF ZACON
  !     IS CALLED FROM ZAIRY.
  !
  !***SEE ALSO  ZAIRY
  !***ROUTINES CALLED  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ZACAI
  !     COMPLEX CSGN,CSPN,C1,C2,Y,Z,ZN,CY
  REAL(8) :: Alim , arg , ascle , az , csgnr , csgni , cspnr , cspni , &
    c1r , c1i , c2r , c2i , cyr , cyi , dfnu , Elim , fmr , &
    Fnu , pi , Rl , sgn , Tol , yy , Yr , Yi , Zr , Zi , &
    znr , zni , D1MACH , ZABS
  INTEGER inu , iuf , Kode , Mr , N , nn , nw , Nz
  DIMENSION Yr(N) , Yi(N) , cyr(2) , cyi(2)
  EXTERNAL ZABS
  DATA pi/3.14159265358979324D0/
  !***FIRST EXECUTABLE STATEMENT  ZACAI
  Nz = 0
  znr = -Zr
  zni = -Zi
  az = ZABS(Zr,Zi)
  nn = N
  dfnu = Fnu + (N-1)
  IF ( az<=2.0D0 ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR THE I FUNCTION
    !-----------------------------------------------------------------------
    CALL ZSERI(znr,zni,Fnu,Kode,nn,Yr,Yi,nw,Tol,Elim,Alim)
  ELSEIF ( az*az*0.25D0>dfnu+1.0D0 ) THEN
    IF ( az<Rl ) THEN
      !-----------------------------------------------------------------------
      !     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL ZMLRI(znr,zni,Fnu,Kode,nn,Yr,Yi,nw,Tol)
      IF ( nw<0 ) GOTO 100
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL ZASYI(znr,zni,Fnu,Kode,nn,Yr,Yi,nw,Rl,Tol,Elim,Alim)
      IF ( nw<0 ) GOTO 100
    ENDIF
  ELSE
    CALL ZSERI(znr,zni,Fnu,Kode,nn,Yr,Yi,nw,Tol,Elim,Alim)
  ENDIF
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
  !-----------------------------------------------------------------------
  CALL ZBKNU(znr,zni,Fnu,Kode,1,cyr,cyi,nw,Tol,Elim,Alim)
  IF ( nw==0 ) THEN
    fmr = Mr
    sgn = -DSIGN(pi,fmr)
    csgnr = 0.0D0
    csgni = sgn
    IF ( Kode/=1 ) THEN
      yy = -zni
      csgnr = -csgni*SIN(yy)
      csgni = csgni*COS(yy)
    ENDIF
    !-----------------------------------------------------------------------
    !     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    !     WHEN FNU IS LARGE
    !-----------------------------------------------------------------------
    inu = Fnu
    arg = (Fnu-inu)*sgn
    cspnr = COS(arg)
    cspni = SIN(arg)
    IF ( MOD(inu,2)/=0 ) THEN
      cspnr = -cspnr
      cspni = -cspni
    ENDIF
    c1r = cyr(1)
    c1i = cyi(1)
    c2r = Yr(1)
    c2i = Yi(1)
    IF ( Kode/=1 ) THEN
      iuf = 0
      ascle = 1.0D+3*D1MACH(1)/Tol
      CALL ZS1S2(znr,zni,c1r,c1i,c2r,c2i,nw,ascle,Alim,iuf)
      Nz = Nz + nw
    ENDIF
    Yr(1) = cspnr*c1r - cspni*c1i + csgnr*c2r - csgni*c2i
    Yi(1) = cspnr*c1i + cspni*c1r + csgnr*c2i + csgni*c2r
    RETURN
  ENDIF
  100  Nz = -1
  IF ( nw==(-2) ) Nz = -2
END SUBROUTINE ZACAI
