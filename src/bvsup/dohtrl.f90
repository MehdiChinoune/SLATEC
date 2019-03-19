!** DOHTRL
SUBROUTINE DOHTRL(Q,N,Nrda,Diag,Irank,Div,Td)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DBVSUP and DSUDS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (OHTROL-S, DOHTRL-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !     For a rank deficient problem, additional orthogonal
  !     HOUSEHOLDER transformations are applied to the left side
  !     of Q to further reduce the triangular form.
  !     Thus, after application of the routines DORTHR and DOHTRL
  !     to the original matrix, the result is a nonsingular
  !     triangular matrix while the remainder of the matrix
  !     has been zeroed out.
  !
  !***
  ! **See also:**  DBVSUP, DSUDS
  !***
  ! **Routines called:**  DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  
  REAL(8) :: DDOT
  INTEGER Irank, irp, j, k, kir, kirm, l, N, nmir, Nrda
  REAL(8) :: dd, Diag(*), diagk, Div(*), Q(Nrda,*), qs, sig, &
    sqd, Td(*), tdv
  !* FIRST EXECUTABLE STATEMENT  DOHTRL
  nmir = N - Irank
  irp = Irank + 1
  DO k = 1, Irank
    kir = irp - k
    diagk = Diag(kir)
    sig = (diagk*diagk) + DDOT(nmir,Q(irp,kir),1,Q(irp,kir),1)
    dd = SIGN(SQRT(sig),-diagk)
    Div(kir) = dd
    tdv = diagk - dd
    Td(kir) = tdv
    IF ( k/=Irank ) THEN
      kirm = kir - 1
      sqd = dd*diagk - sig
      DO j = 1, kirm
        qs = ((tdv*Q(kir,j))+DDOT(nmir,Q(irp,j),1,Q(irp,kir),1))/sqd
        Q(kir,j) = Q(kir,j) + qs*tdv
        DO l = irp, N
          Q(l,j) = Q(l,j) + qs*Q(l,kir)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
END SUBROUTINE DOHTRL
