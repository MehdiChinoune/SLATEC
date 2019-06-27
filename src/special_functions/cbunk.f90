!** CBUNK
PURE SUBROUTINE CBUNK(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  !> Subsidiary to CBESH and CBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBUNK-A, ZBUNK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU>FNUL.
  !     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
  !     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2
  !
  !***
  ! **See also:**  CBESH, CBESK
  !***
  ! **Routines called:**  CUNK1, CUNK2

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(SP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(SP), INTENT(IN) :: Z
  COMPLEX(SP), INTENT(OUT) :: Y(N)
  REAL(SP) :: ax, ay, xx, yy
  !* FIRST EXECUTABLE STATEMENT  CBUNK
  Nz = 0
  xx = REAL(Z)
  yy = AIMAG(Z)
  ax = ABS(xx)*1.7321_SP
  ay = ABS(yy)
  IF( ay>ax ) THEN
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
    !     APPLIED IN PI/3<ABS(ARG(Z))<=PI/2 WHERE M=+I OR -I
    !     AND HPI=PI/2
    !-----------------------------------------------------------------------
    CALL CUNK2(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  ELSE
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
    !     -PI/3<=ARG(Z)<=PI/3
    !-----------------------------------------------------------------------
    CALL CUNK1(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  END IF
END SUBROUTINE CBUNK