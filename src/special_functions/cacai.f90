!** CACAI
PURE SUBROUTINE CACAI(Z,Fnu,Kode,Mr,N,Y,Nz,Rl,Tol,Elim,Alim)
  !> Subsidiary to CAIRY
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CACAI-A, ZACAI-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
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
  !***
  ! **See also:**  CAIRY
  !***
  ! **Routines called:**  CASYI, CBKNU, CMLRI, CS1S2, CSERI, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_sp
  !
  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alim, Elim, Fnu, Rl, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Y(N)
  !
  INTEGER :: inu, iuf, nn, nw
  COMPLEX(SP) :: csgn, cspn, c1, c2, zn, cy(2)
  REAL(SP) :: arg, ascle, az, cpn, dfnu, fmr, sgn, spn, yy
  REAL(SP), PARAMETER :: pi = 3.14159265358979324_SP
  !* FIRST EXECUTABLE STATEMENT  CACAI
  Nz = 0
  zn = -Z
  az = ABS(Z)
  nn = N
  dfnu = Fnu + (N-1)
  IF( az<=2._SP ) THEN
    !-----------------------------------------------------------------------
    !     POWER SERIES FOR THE I FUNCTION
    !-----------------------------------------------------------------------
    CALL CSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  ELSEIF( az*az*0.25E0>dfnu+1._SP ) THEN
    IF( az<Rl ) THEN
      !-----------------------------------------------------------------------
      !     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL CMLRI(zn,Fnu,Kode,nn,Y,nw,Tol)
      IF( nw<0 ) GOTO 100
    ELSE
      !-----------------------------------------------------------------------
      !     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
      !-----------------------------------------------------------------------
      CALL CASYI(zn,Fnu,Kode,nn,Y,nw,Rl,Tol,Elim,Alim)
      IF( nw<0 ) GOTO 100
    END IF
  ELSE
    CALL CSERI(zn,Fnu,Kode,nn,Y,nw,Tol,Elim,Alim)
  END IF
  !-----------------------------------------------------------------------
  !     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
  !-----------------------------------------------------------------------
  CALL CBKNU(zn,Fnu,Kode,1,cy,nw,Tol,Elim,Alim)
  IF( nw==0 ) THEN
    fmr = Mr
    sgn = -SIGN(pi,fmr)
    csgn = CMPLX(0._SP,sgn,SP)
    IF( Kode/=1 ) THEN
      yy = -AIMAG(zn)
      cpn = COS(yy)
      spn = SIN(yy)
      csgn = csgn*CMPLX(cpn,spn,SP)
    END IF
    !-----------------------------------------------------------------------
    !     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
    !     WHEN FNU IS LARGE
    !-----------------------------------------------------------------------
    inu = INT( Fnu )
    arg = (Fnu-inu)*sgn
    cpn = COS(arg)
    spn = SIN(arg)
    cspn = CMPLX(cpn,spn,SP)
    IF( MOD(inu,2)==1 ) cspn = -cspn
    c1 = cy(1)
    c2 = Y(1)
    IF( Kode/=1 ) THEN
      iuf = 0
      ascle = 1.E+3_SP*tiny_sp/Tol
      CALL CS1S2(zn,c1,c2,nw,ascle,Alim,iuf)
      Nz = Nz + nw
    END IF
    Y(1) = cspn*c1 + csgn*c2
    RETURN
  END IF
  100  Nz = -1
  IF( nw==(-2) ) Nz = -2
  !
END SUBROUTINE CACAI