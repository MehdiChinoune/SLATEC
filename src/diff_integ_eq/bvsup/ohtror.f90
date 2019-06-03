!** OHTROR
SUBROUTINE OHTROR(Q,N,Nrda,Diag,Irank,Div,Td)
  !>
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (OHTROR-S)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !     For a rank deficient problem, additional orthogonal
  !     HOUSEHOLDER transformations are applied to the right side
  !     of Q to further reduce the triangular form.
  !     Thus, after application of the routines ORTHOL and OHTROR
  !     to the original matrix, the result is a nonsingular
  !     triangular matrix while the remainder of the matrix
  !     has been zeroed out.
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  SDOT

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  INTEGER :: Irank, N, Nrda
  REAL(SP) :: Diag(Irank), Div(Irank), Q(Nrda,Irank), Td(Irank)
  INTEGER :: irp, j, k, kir, kirm, l, nmir
  REAL(SP) :: dd, diagk, qs, sig, sqd, tdv
  !* FIRST EXECUTABLE STATEMENT  OHTROR
  nmir = N - Irank
  irp = Irank + 1
  DO k = 1, Irank
    kir = irp - k
    diagk = Diag(kir)
    sig = diagk**2 + NORM2(Q(kir,irp:N))**2
    dd = SIGN(SQRT(sig),-diagk)
    Div(kir) = dd
    tdv = diagk - dd
    Td(kir) = tdv
    IF ( k/=Irank ) THEN
      kirm = kir - 1
      sqd = dd*diagk - sig
      DO j = 1, kirm
        qs = ((tdv*Q(j,kir))+DOT_PRODUCT(Q(j,irp:N),Q(kir,irp:N)))/sqd
        Q(j,kir) = Q(j,kir) + qs*tdv
        DO l = irp, N
          Q(j,l) = Q(j,l) + qs*Q(kir,l)
        END DO
      END DO
    END IF
  END DO
END SUBROUTINE OHTROR
