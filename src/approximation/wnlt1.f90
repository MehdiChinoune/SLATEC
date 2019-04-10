!** WNLT1
SUBROUTINE WNLT1(I,Lend,Mend,Ir,Mdw,Recalc,Imax,Hbar,H,Scalee,W)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to WNLIT
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (WNLT1-S, DWNLT1-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***
  ! **Description:**
  !
  !     To update the column Sum Of Squares and find the pivot column.
  !     The column Sum of Squares Vector will be updated at each step.
  !     When numerically necessary, these values will be recomputed.
  !
  !***
  ! **See also:**  WNLIT
  !***
  ! **Routines called:**  ISAMAX

  !* REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLIT and made a subroutine.  (RWC))

  INTEGER I, Imax, Ir, Lend, Mdw, Mend
  REAL H(*), Hbar, Scalee(*), W(Mdw,*)
  LOGICAL Recalc
  !
  INTEGER, EXTERNAL :: ISAMAX
  !
  INTEGER j, k
  !
  !* FIRST EXECUTABLE STATEMENT  WNLT1
  IF ( Ir/=1.AND.(.NOT.Recalc) ) THEN
    !
    !        Update column SS=summ of squares.
    !
    DO j = I, Lend
      H(j) = H(j) - Scalee(Ir-1)*W(Ir-1,j)**2
    END DO
    !
    !        Test for numerical accuracy.
    !
    Imax = ISAMAX(Lend-I+1,H(I),1) + I - 1
    Recalc = (Hbar+1.E-3*H(Imax))==Hbar
  END IF
  !
  !     If required, recalculate column SS, using rows IR through MEND.
  !
  IF ( Recalc ) THEN
    DO j = I, Lend
      H(j) = 0.E0
      DO k = Ir, Mend
        H(j) = H(j) + Scalee(k)*W(k,j)**2
      END DO
    END DO
    !
    !        Find column with largest SS.
    !
    Imax = ISAMAX(Lend-I+1,H(I),1) + I - 1
    Hbar = H(Imax)
  END IF
END SUBROUTINE WNLT1
