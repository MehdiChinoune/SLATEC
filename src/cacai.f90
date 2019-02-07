!*==CACAI.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CACAI
SUBROUTINE CACAI(Z,Fnu,Kode,Mr,N,Y,Nz,Rl,Tol,Elim,Alim)
  IMPLICIT NONE
  !*--CACAI5
  !***BEGIN PROLOGUE  CACAI
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CAIRY
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CACAI-A, ZACAI-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA
  !
  !         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
  !                 MP=PI*MR*CMPLX(0.0,1.0)
  !
  !     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
  !     HALF Z PLANE FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
  !     CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
  !     RECURRENCE REMOVED. A RECURSIVE CALL TO CACON CAN RESULT IF CACON
  !     IS CALLED FROM CAIRY.
  !
  !***SEE ALSO  CAIRY
  !***ROUTINES CALLED  CASYI, CBKNU, CMLRI, CS1S2, CSERI, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CACAI
  COMPLEX csgn, cspn, c1, c2, Y, Z, zn, cy
  REAL Alim, arg, ascle, az, cpn, dfnu, Elim, fmr, Fnu, pi, Rl, &
    sgn, spn, Tol, yy, R1MACH
  INTEGER inu, iuf, Kode, Mr, N, nn, nw, Nz
  DIMENSION Y(N), cy(2)
  DATA pi/3.14159265358979324E0/
  !***FIRST EXECUTABLE STATEMENT  CACAI
  Nz = 0
  zn = -Z
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)
  IF ( az<=2.0E0 ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR THE I FUNCTION
    !-----------------------------------------------------------------------
    CALL CSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  ELSEIF ( az*az*0.25E0>dfnu+1.0E0 ) THEN
    IF ( az<Rl ) THEN
      !-----------------------------------------------------------------------
      !     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL CMLRI(zn,Fnu,Kode,nn,Y,nw,Tol)
      IF ( nw<0 ) GOTO 100
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL CASYI(zn,Fnu,Kode,nn,Y,nw,Rl,Tol,Elim,Alim)
      IF ( nw<0 ) GOTO 100
    ENDIF
  ELSE
    CALL CSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  ENDIF
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
  !-----------------------------------------------------------------------
  CALL CBKNU(zn,Fnu,Kode,1,cy,nw,Tol,Elim,Alim)
  IF ( nw==0 ) THEN
    fmr = Mr
    sgn = -SIGN(pi,fmr)
    csgn = CMPLX(0.0E0,sgn)
    IF ( Kode/=1 ) THEN
      yy = -AIMAG(zn)
      cpn = COS(yy)
      spn = SIN(yy)
      csgn = csgn*CMPLX(cpn,spn)
    ENDIF
    !-----------------------------------------------------------------------
    !     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    !     WHEN FNU IS LARGE
    !-----------------------------------------------------------------------
    inu = Fnu
    arg = (Fnu-inu)*sgn
    cpn = COS(arg)
    spn = SIN(arg)
    cspn = CMPLX(cpn,spn)
    IF ( MOD(inu,2)==1 ) cspn = -cspn
    c1 = cy(1)
    c2 = Y(1)
    IF ( Kode/=1 ) THEN
      iuf = 0
      ascle = 1.0E+3*R1MACH(1)/Tol
      CALL CS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
      Nz = Nz + nw
    ENDIF
    Y(1) = cspn*c1 + csgn*c2
    RETURN
  ENDIF
  100  Nz = -1
  IF ( nw==(-2) ) Nz = -2
END SUBROUTINE CACAI
