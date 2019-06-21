!** ORTHO4
SUBROUTINE ORTHO4(Usol,Idmn,Zn,Zm,Pertrb)
  !> Subsidiary to SEPX4
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
  USE SPL4, ONLY : is_com, js_com, ms_com, ns_com
  INTEGER :: Idmn
  REAL(SP) :: Pertrb, Usol(Idmn,ns_com), Zm(ms_com), Zn(ns_com)
  INTEGER :: i, ifnl, ii, istr, j, jfnl, jj, jstr
  REAL(SP) :: ete, ute
  !* FIRST EXECUTABLE STATEMENT  ORTHO4
  istr = is_com
  ifnl = ms_com
  jstr = js_com
  jfnl = ns_com
  !
  !     COMPUTE WEIGHTED INNER PRODUCTS
  !
  ute = 0._SP
  ete = 0._SP
  DO i = is_com, ms_com
    ii = i - is_com + 1
    DO j = js_com, ns_com
      jj = j - js_com + 1
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
