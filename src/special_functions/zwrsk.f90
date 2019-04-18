!** ZWRSK
SUBROUTINE ZWRSK(Zrr,Zri,Fnu,Kode,N,Yr,Yi,Nz,Cwr,Cwi,Tol,Elim,Alim)
  !>
  !***
  !  Subsidiary to ZBESI and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CWRSK-A, ZWRSK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
  !     NORMALIZING THE I FUNCTION RATIOS FROM ZRATI BY THE WRONSKIAN
  !
  !***
  ! **See also:**  ZBESI, ZBESK
  !***
  ! **Routines called:**  D1MACH, ZABS, ZBKNU, ZRATI

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  !     COMPLEX CINU,CSCL,CT,CW,C1,C2,RCT,ST,Y,ZR
  INTEGER i, Kode, N, nw, Nz
  REAL(8) :: act, acw, Alim, ascle, cinui, cinur, csclr, cti, &
    ctr, Cwi(2), Cwr(2), c1i, c1r, c2i, c2r, Elim, Fnu, &
    pti, ptr, ract, sti, str, Tol, Yi(N), Yr(N), Zri, Zrr
  !* FIRST EXECUTABLE STATEMENT  ZWRSK
  !-----------------------------------------------------------------------
  !     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
  !     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
  !     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
  !-----------------------------------------------------------------------
  !
  Nz = 0
  CALL ZBKNU(Zrr,Zri,Fnu,Kode,2,Cwr,Cwi,nw,Tol,Elim,Alim)
  IF ( nw/=0 ) THEN
    Nz = -1
    IF ( nw==(-2) ) Nz = -2
    RETURN
  ELSE
    CALL ZRATI(Zrr,Zri,Fnu,N,Yr,Yi,Tol)
    !-----------------------------------------------------------------------
    !     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    !     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    !-----------------------------------------------------------------------
    cinur = 1.0D0
    cinui = 0.0D0
    IF ( Kode/=1 ) THEN
      cinur = COS(Zri)
      cinui = SIN(Zri)
    END IF
    !-----------------------------------------------------------------------
    !     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    !     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    !     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    !     THE RESULT IS ON SCALE.
    !-----------------------------------------------------------------------
    acw = ZABS(Cwr(2),Cwi(2))
    ascle = 1.0D+3*D1MACH(1)/Tol
    csclr = 1.0D0
    IF ( acw>ascle ) THEN
      ascle = 1.0D0/ascle
      IF ( acw>=ascle ) csclr = Tol
    ELSE
      csclr = 1.0D0/Tol
    END IF
  END IF
  c1r = Cwr(1)*csclr
  c1i = Cwi(1)*csclr
  c2r = Cwr(2)*csclr
  c2i = Cwi(2)*csclr
  str = Yr(1)
  sti = Yi(1)
  !-----------------------------------------------------------------------
  !     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0D0/ABS(CT) PREVENTS
  !     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
  !-----------------------------------------------------------------------
  ptr = str*c1r - sti*c1i
  pti = str*c1i + sti*c1r
  ptr = ptr + c2r
  pti = pti + c2i
  ctr = Zrr*ptr - Zri*pti
  cti = Zrr*pti + Zri*ptr
  act = ZABS(ctr,cti)
  ract = 1.0D0/act
  ctr = ctr*ract
  cti = -cti*ract
  ptr = cinur*ract
  pti = cinui*ract
  cinur = ptr*ctr - pti*cti
  cinui = ptr*cti + pti*ctr
  Yr(1) = cinur*csclr
  Yi(1) = cinui*csclr
  IF ( N==1 ) RETURN
  DO i = 2, N
    ptr = str*cinur - sti*cinui
    cinui = str*cinui + sti*cinur
    cinur = ptr
    str = Yr(i)
    sti = Yi(i)
    Yr(i) = cinur*csclr
    Yi(i) = cinui*csclr
  END DO
  RETURN
END SUBROUTINE ZWRSK
