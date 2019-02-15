!DECK WNLT1
SUBROUTINE WNLT1(I,Lend,Mend,Ir,Mdw,Recalc,Imax,Hbar,H,Scale,W)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  WNLT1
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to WNLIT
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (WNLT1-S, DWNLT1-D)
  !***AUTHOR  Hanson, R. J., (SNLA)
  !           Haskell, K. H., (SNLA)
  !***DESCRIPTION
  !
  !     To update the column Sum Of Squares and find the pivot column.
  !     The column Sum of Squares Vector will be updated at each step.
  !     When numerically necessary, these values will be recomputed.
  !
  !***SEE ALSO  WNLIT
  !***ROUTINES CALLED  ISAMAX
  !***REVISION HISTORY  (YYMMDD)
  !   790701  DATE WRITTEN
  !   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
  !***END PROLOGUE  WNLT1
  INTEGER I, Imax, Ir, Lend, Mdw, Mend
  REAL H(*), Hbar, Scale(*), W(Mdw,*)
  LOGICAL Recalc
  !
  EXTERNAL ISAMAX
  INTEGER ISAMAX
  !
  INTEGER j, k
  !
  !***FIRST EXECUTABLE STATEMENT  WNLT1
  IF ( Ir/=1.AND.(.NOT.Recalc) ) THEN
    !
    !        Update column SS=sum of squares.
    !
    DO j = I, Lend
      H(j) = H(j) - Scale(Ir-1)*W(Ir-1,j)**2
    ENDDO
    !
    !        Test for numerical accuracy.
    !
    Imax = ISAMAX(Lend-I+1,H(I),1) + I - 1
    Recalc = (Hbar+1.E-3*H(Imax))==Hbar
  ENDIF
  !
  !     If required, recalculate column SS, using rows IR through MEND.
  !
  IF ( Recalc ) THEN
    DO j = I, Lend
      H(j) = 0.E0
      DO k = Ir, Mend
        H(j) = H(j) + Scale(k)*W(k,j)**2
      ENDDO
    ENDDO
    !
    !        Find column with largest SS.
    !
    Imax = ISAMAX(Lend-I+1,H(I),1) + I - 1
    Hbar = H(Imax)
  ENDIF
END SUBROUTINE WNLT1
