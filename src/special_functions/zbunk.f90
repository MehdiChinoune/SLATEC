!** ZBUNK
SUBROUTINE ZBUNK(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Tol,Elim,Alim)
  !> Subsidiary to ZBESH and ZBESK
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (CBUNI-A, ZBUNI-A)
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

  !     COMPLEX Y,Z
  INTEGER :: Kode, Mr, N, Nz
  REAL(DP) :: Alim, ax, ay, Elim, Fnu, Tol, Yi(N), Yr(N), Zi, Zr
  !* FIRST EXECUTABLE STATEMENT  ZBUNK
  Nz = 0
  ax = ABS(Zr)*1.7321_DP
  ay = ABS(Zi)
  IF( ay>ax ) THEN
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
    !     APPLIED IN PI/3<ABS(ARG(Z))<=PI/2 WHERE M=+I OR -I
    !     AND HPI=PI/2
    !-----------------------------------------------------------------------
    CALL ZUNK2(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Tol,Elim,Alim)
  ELSE
    !-----------------------------------------------------------------------
    !     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
    !     -PI/3<=ARG(Z)<=PI/3
    !-----------------------------------------------------------------------
    CALL ZUNK1(Zr,Zi,Fnu,Kode,Mr,N,Yr,Yi,Nz,Tol,Elim,Alim)
  END IF
END SUBROUTINE ZBUNK
