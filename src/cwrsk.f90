!DECK CWRSK
SUBROUTINE CWRSK(Zr,Fnu,Kode,N,Y,Nz,Cw,Tol,Elim,Alim)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CWRSK
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBESI and CBESK
  !***LIBRARY   SLATEC
  !***TYPE      ALL (CWRSK-A, ZWRSK-A)
  !***AUTHOR  Amos, D. E., (SNL)
  !***DESCRIPTION
  !
  !     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
  !     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN
  !
  !***SEE ALSO  CBESI, CBESK
  !***ROUTINES CALLED  CBKNU, CRATI, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CWRSK
  COMPLEX cinu, cscl, ct, Cw, c1, c2, rct, st, Y, Zr
  REAL act, acw, Alim, ascle, Elim, Fnu, s1, s2, Tol, yy, R1MACH
  INTEGER i, Kode, N, nw, Nz
  DIMENSION Y(N), Cw(2)
  !***FIRST EXECUTABLE STATEMENT  CWRSK
  !-----------------------------------------------------------------------
  !     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
  !     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
  !     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
  !-----------------------------------------------------------------------
  Nz = 0
  CALL CBKNU(Zr,Fnu,Kode,2,Cw,nw,Tol,Elim,Alim)
  IF ( nw/=0 ) THEN
    Nz = -1
    IF ( nw==(-2) ) Nz = -2
    GOTO 99999
  ELSE
    CALL CRATI(Zr,Fnu,N,Y,Tol)
    !-----------------------------------------------------------------------
    !     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
    !     R(FNU+J-1,Z)=Y(J),  J=1,...,N
    !-----------------------------------------------------------------------
    cinu = CMPLX(1.0E0,0.0E0)
    IF ( Kode/=1 ) THEN
      yy = AIMAG(Zr)
      s1 = COS(yy)
      s2 = SIN(yy)
      cinu = CMPLX(s1,s2)
    ENDIF
    !-----------------------------------------------------------------------
    !     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
    !     THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
    !     SCALED TO PREVENT OVER OR UNDERFLOW. CUOIK HAS DETERMINED THAT
    !     THE RESULT IS ON SCALE.
    !-----------------------------------------------------------------------
    acw = ABS(Cw(2))
    ascle = 1.0E+3*R1MACH(1)/Tol
    cscl = CMPLX(1.0E0,0.0E0)
    IF ( acw>ascle ) THEN
      ascle = 1.0E0/ascle
      IF ( acw>=ascle ) cscl = CMPLX(Tol,0.0E0)
    ELSE
      cscl = CMPLX(1.0E0/Tol,0.0E0)
    ENDIF
  ENDIF
  c1 = Cw(1)*cscl
  c2 = Cw(2)*cscl
  st = Y(1)
  !-----------------------------------------------------------------------
  !     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0E0/ABS(CT) PREVENTS
  !     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
  !-----------------------------------------------------------------------
  ct = Zr*(c2+st*c1)
  act = ABS(ct)
  rct = CMPLX(1.0E0/act,0.0E0)
  ct = CONJG(ct)*rct
  cinu = cinu*rct*ct
  Y(1) = cinu*cscl
  IF ( N==1 ) RETURN
  DO i = 2, N
    cinu = st*cinu
    st = Y(i)
    Y(i) = cinu*cscl
  ENDDO
  RETURN
  99999 CONTINUE
  END SUBROUTINE CWRSK
