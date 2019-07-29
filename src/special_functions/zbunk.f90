!** ZBUNK
PURE SUBROUTINE ZBUNK(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZBESH and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBUNK-A, ZBUNK-A)
  !***
  ! **Author:**  Amos, D. E., (SNL)
  !***
  ! **Description:**
  !
  !     ZBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU>FNUL.
  !     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
  !     IN ZUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN ZUNK2
  !
  !***
  ! **See also:**  ZBESH, ZBESK
  !***
  ! **Routines called:**  ZUNK1, ZUNK2

  !* REVISION HISTORY  (YYMMDD)
  !   830501  DATE WRITTEN
  !   910415  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER, INTENT(IN) :: Kode, Mr, N
  INTEGER, INTENT(OUT) :: Nz
  REAL(DP), INTENT(IN) :: Alim, Elim, Fnu, Tol
  COMPLEX(DP), INTENT(IN) :: Z
  COMPLEX(DP), INTENT(OUT) :: Y(N)
  !
  REAL(DP) :: ax, ay, xx, yy
  !* FIRST EXECUTABLE STATEMENT  ZBUNK
  Nz = 0
  xx = REAL(Z,DP)
  yy = AIMAG(Z)
  ax = ABS(xx)*1.7321_DP
  ay = ABS(yy)
  IF( ay>ax ) THEN
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
    !     APPLIED IN PI/3<ABS(ARG(Z))<=PI/2 WHERE M=+I OR -I AND HPI=PI/2
    !-----------------------------------------------------------------------
    CALL ZUNK2(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  ELSE
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
    !     -PI/3<=ARG(Z)<=PI/3
    !-----------------------------------------------------------------------
    CALL ZUNK1(Z,Fnu,Kode,Mr,N,Y,Nz,Tol,Elim,Alim)
  END IF
  !
END SUBROUTINE ZBUNK