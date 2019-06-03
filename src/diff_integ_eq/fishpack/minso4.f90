!** MINSO4
SUBROUTINE MINSO4(Usol,Idmn,Zn,Zm)
  !>
  !  Subsidiary to SEPX4
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (MINSO4-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine orthogonalizes the array USOL with respect to
  !     the constant array in a weighted least squares norm.
  !
  !     Entry at MINSO4 occurs when the final solution is
  !     to be minimized with respect to the weighted
  !     least squares norm.
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
  USE SPL4, ONLY : l_com, is_com, js_com, k_com, ms_com, ns_com
  INTEGER :: Idmn
  REAL(SP) :: Usol(Idmn,l_com), Zm(ms_com), Zn(ns_com)
  INTEGER :: i, ifnl, ii, istr, j, jfnl, jj, jstr
  REAL(SP) :: ete, pertrb, ute
  !* FIRST EXECUTABLE STATEMENT  MINSO4
  istr = 1
  ifnl = k_com
  jstr = 1
  jfnl = l_com
  !
  !     COMPUTE WEIGHTED INNER PRODUCTS
  !
  ute = 0.0
  ete = 0.0
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
  pertrb = ute/ete
  !
  !     SUBTRACT OFF CONSTANT PERTRB
  !
  DO i = istr, ifnl
    DO j = jstr, jfnl
      Usol(i,j) = Usol(i,j) - pertrb
    END DO
  END DO
END SUBROUTINE MINSO4
