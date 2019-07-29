!** ZWRSK
PURE SUBROUTINE ZWRSK(Zr,Fnu,Kode,N,Y,Nz,Cw,Tol,Elim,Alim)
  !> Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CWRSK-A, ZWRSK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z)>=0.0 BY
  !     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZBKNU, ZRATI

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : tiny_dp
  !
  INTEGER, INTENT(IN) :: Kode, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Zr
  COMPLEX(DP), INTENT(OUT) :: Cw(2), Y(N)
  !
  INTEGER :: i, nw
  COMPLEX(DP) :: cinu, cscl, ct, c1, c2, rct, st
  REAL(DP) :: act, acw, ascle, s1, s2, yy
  !* FIRST EXECUTABLE STATEMENT  ZWRSK
  !-----------------------------------------------------------------------
  !     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
  !     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM ZRATI NORMALIZED BY THE
  !     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM ZBKNU.
  !-----------------------------------------------------------------------
  Nz = 0
  CALL ZBKNU(Zr,Fnu,Kode,2,Cw,nw,Tol,Elim,Alim)
  IF( nw/=0 ) THEN
    Nz = -1
    IF( nw==(-2) ) Nz = -2
    RETURN
  ELSE
    CALL ZRATI(Zr,Fnu,N,Y,Tol)
    !-----------------------------------------------------------------------
    !     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    !     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    !-----------------------------------------------------------------------
    cinu = CMPLX(1._DP,0._DP,DP)
    IF( Kode/=1 ) THEN
      yy = AIMAG(Zr)
      s1 = COS(yy)
      s2 = SIN(yy)
      cinu = CMPLX(s1,s2,DP)
    END IF
    !-----------------------------------------------------------------------
    !     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    !     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    !     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    !     THE RESULT IS ON SCALE.
    !-----------------------------------------------------------------------
    acw = ABS(Cw(2))
    ascle = 1.E+3_DP*tiny_dp/Tol
    cscl = CMPLX(1._DP,0._DP,DP)
    IF( acw>ascle ) THEN
      ascle = 1._DP/ascle
      IF( acw>=ascle ) cscl = CMPLX(Tol,0._DP,DP)
    ELSE
      cscl = CMPLX(1._DP/Tol,0._DP,DP)
    END IF
  END IF
  c1 = Cw(1)*cscl
  c2 = Cw(2)*cscl
  st = Y(1)
  !-----------------------------------------------------------------------
  !     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0E0/ABS(CT) PREVENTS
  !     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
  !-----------------------------------------------------------------------
  ct = Zr*(c2+st*c1)
  act = ABS(ct)
  rct = CMPLX(1._DP/act,0._DP,DP)
  ct = CONJG(ct)*rct
  cinu = cinu*rct*ct
  Y(1) = cinu*cscl
  IF( N==1 ) RETURN
  DO i = 2, N
    cinu = st*cinu
    st = Y(i)
    Y(i) = cinu*cscl
  END DO
  !
  RETURN
END SUBROUTINE ZWRSK