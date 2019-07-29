!** ZACAI
PURE SUBROUTINE ZACAI(Z,Fnu,Kode,Mr,N,Y,Nz,Rl,Tol,Elim,Alim)
  !> Subsidiary to ZAIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CACAI-A, ZACAI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
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
  !***
  ! **See also:**  ZAIRY
  !***
  ! **Routines called:**  D1MACH, ZABS, ZASYI, ZBKNU, ZMLRI, ZS1S2, ZSERI

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_dp
  !
  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Rl, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: inu, iuf, nn, nw
  COMPLEX(DP) :: csgn, cspn, c1, c2, zn, cy(2)
  REAL(DP) :: arg, ascle, az, cpn, dfnu, fmr, sgn, spn, yy
  REAL(DP), PARAMETER :: pi = 3.14159265358979324_DP
  !* FIRST EXECUTABLE STATEMENT  ZACAI
  Nz = 0
  zn = -Z
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)
  IF( az<=2._DP ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR THE I FUNCTION
    !-----------------------------------------------------------------------
    CALL ZSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  ELSEIF( az*az*0.25_DP>dfnu+1._DP ) THEN
    IF( az<Rl ) THEN
      !-----------------------------------------------------------------------
      !     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL ZMLRI(zn,Fnu,Kode,nn,Y,nw,Tol)
      IF( nw<0 ) GOTO 100
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL ZASYI(zn,Fnu,Kode,nn,Y,nw,Rl,Tol,Elim,Alim)
      IF( nw<0 ) GOTO 100
    END IF
  ELSE
    CALL ZSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  END IF
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
  !-----------------------------------------------------------------------
  CALL ZBKNU(zn,Fnu,Kode,1,cy,nw,Tol,Elim,Alim)
  IF( nw==0 ) THEN
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
    c1 = cy(1)
    c2 = Y(1)
    IF( Kode/=1 ) THEN
      iuf = 0
      ascle = 1.E+3_DP*tiny_dp/Tol
      CALL ZS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
      Nz = Nz + nw
    END IF
    Y(1) = cspn*c1 + csgn*c2
    RETURN
  END IF
  100  Nz = -1
  IF( nw==(-2) ) Nz = -2
  !
END SUBROUTINE ZACAI