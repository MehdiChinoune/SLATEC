!** ORTHO4
SUBROUTINE ORTHO4(Usol,Idmn,Zn,Zm,Pertrb)
  !>
  !  Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (ORTHO4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine orthogonalizes the array USOL with respect to
  !     the constant array in a weighted least squares norm.
  !
  !***
  ! **See also:**  SEPX4
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPL4

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPL4, ONLY : IS, JS, MS, NS
  INTEGER i, Idmn, ifnl, ii, istr, j, jfnl, jj, jstr
  REAL ete, Pertrb, Usol(Idmn,*), ute, Zm(*), Zn(*)
  !* FIRST EXECUTABLE STATEMENT  ORTHO4
  istr = IS
  ifnl = MS
  jstr = JS
  jfnl = NS
  !
  !     COMPUTE WEIGHTED INNER PRODUCTS
  !
  ute = 0.0
  ete = 0.0
  DO i = IS, MS
    ii = i - IS + 1
    DO j = JS, NS
      jj = j - JS + 1
      ete = ete + Zm(ii)*Zn(jj)
      ute = ute + Usol(i,j)*Zm(ii)*Zn(jj)
    END DO
  END DO
  !
  !     SET PERTURBATION PARAMETER
  !
  Pertrb = ute/ete
  !
  !     SUBTRACT OFF CONSTANT PERTRB
  !
  DO i = istr, ifnl
    DO j = jstr, jfnl
      Usol(i,j) = Usol(i,j) - Pertrb
    END DO
  END DO
END SUBROUTINE ORTHO4
