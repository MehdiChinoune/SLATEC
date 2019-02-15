!DECK OHTROR
SUBROUTINE OHTROR(Q,N,Nrda,Diag,Irank,Div,Td)
  IMPLICIT NONE
  REAL dd, Diag, diagk, Div, Q, qs, SDOT, sig, sqd, Td, tdv
  INTEGER Irank, irp, j, k, kir, kirm, l, N, nmir, Nrda
  !***BEGIN PROLOGUE  OHTROR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (OHTROR-S)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     For a rank deficient problem, additional orthogonal
  !     HOUSEHOLDER transformations are applied to the right side
  !     of Q to further reduce the triangular form.
  !     Thus, after application of the routines ORTHOL and OHTROR
  !     to the original matrix, the result is a nonsingular
  !     triangular matrix while the remainder of the matrix
  !     has been zeroed out.
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  SDOT
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  OHTROR
  DIMENSION Q(Nrda,*), Diag(*), Div(*), Td(*)
  !***FIRST EXECUTABLE STATEMENT  OHTROR
  nmir = N - Irank
  irp = Irank + 1
  DO k = 1, Irank
    kir = irp - k
    diagk = Diag(kir)
    sig = (diagk*diagk) + SDOT(nmir,Q(kir,irp),Nrda,Q(kir,irp),Nrda)
    dd = SIGN(SQRT(sig),-diagk)
    Div(kir) = dd
    tdv = diagk - dd
    Td(kir) = tdv
    IF ( k/=Irank ) THEN
      kirm = kir - 1
      sqd = dd*diagk - sig
      DO j = 1, kirm
        qs = ((tdv*Q(j,kir))+SDOT(nmir,Q(j,irp),Nrda,Q(kir,irp),Nrda))/sqd
        Q(j,kir) = Q(j,kir) + qs*tdv
        DO l = irp, N
          Q(j,l) = Q(j,l) + qs*Q(kir,l)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
END SUBROUTINE OHTROR
