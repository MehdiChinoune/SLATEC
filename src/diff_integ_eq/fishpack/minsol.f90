!** MINSOL
SUBROUTINE MINSOL(Usol,Idmn,Zn,Zm)
  !>
  !  Subsidiary to SEPELI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (MINSOL-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine orthogonalizes the array USOL with respect to
  !     the constant array in a weighted least squares norm.
  !
  !     Entry at MINSOL occurs when the final solution is
  !     to be minimized with respect to the weighted
  !     least squares norm.
  !
  !***
  ! **See also:**  SEPELI
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    SPLPCM

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  USE SPLPCM, ONLY : L, IS, JS, K, MS, NS
  INTEGER :: Idmn
  REAL :: Usol(Idmn,L), Zm(MS), Zn(NS)
  INTEGER :: i, ifnl, ii, istr, j, jfnl, jj, jstr
  REAL :: ete, pertrb, ute
  !* FIRST EXECUTABLE STATEMENT  MINSOL
  istr = 1
  ifnl = K
  jstr = 1
  jfnl = L
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
  pertrb = ute/ete
  !
  !     SUBTRACT OFF CONSTANT PERTRB
  !
  DO i = istr, ifnl
    DO j = jstr, jfnl
      Usol(i,j) = Usol(i,j) - pertrb
    END DO
  END DO
END SUBROUTINE MINSOL
