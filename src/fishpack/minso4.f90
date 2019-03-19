!** MINSO4
SUBROUTINE MINSO4(Usol,Idmn,Zn,Zm,Pertb)
  IMPLICIT NONE
  !>
  !***
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
  
  REAL AIT, BIT, CIT, DIT, DLX, DLX4, DLY, DLY4, ete, Pertb, &
    pertrb, TDLx3, TDLy3, Usol, ute, Zm, Zn
  INTEGER i, Idmn, ifnl, ii, IS, istr, j, jfnl, jj, JS, jstr, K, &
    KSWx, KSWy, L, MIT, MS, NIT, NS
  COMMON /SPL4  / KSWx, KSWy, K, L, AIT, BIT, CIT, DIT, MIT, NIT, &
    IS, MS, JS, NS, DLX, DLY, TDLx3, TDLy3, DLX4, DLY4
  DIMENSION Usol(Idmn,*), Zn(*), Zm(*)
  !* FIRST EXECUTABLE STATEMENT  MINSO4
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
    ENDDO
  ENDDO
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
    ENDDO
  ENDDO
END SUBROUTINE MINSO4
