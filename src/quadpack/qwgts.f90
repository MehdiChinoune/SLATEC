!DECK QWGTS
REAL FUNCTION QWGTS(X,A,B,Alfa,Beta,Integr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QWGTS
  !***SUBSIDIARY
  !***PURPOSE  This function subprogram is used together with the
  !            routine QAWS and defines the WEIGHT function.
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QWGTS-S, DQWGTS-D)
  !***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
  !             WEIGHT FUNCTION
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***SEE ALSO  QK15W
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  QWGTS
  !
  REAL A, Alfa, B, Beta, bmx, X, xma
  INTEGER Integr
  !***FIRST EXECUTABLE STATEMENT  QWGTS
  xma = X - A
  bmx = B - X
  QWGTS = xma**Alfa*bmx**Beta
  SELECT CASE (Integr)
    CASE (1)
    CASE (3)
      QWGTS = QWGTS*LOG(bmx)
    CASE (4)
      QWGTS = QWGTS*LOG(xma)*LOG(bmx)
    CASE DEFAULT
      QWGTS = QWGTS*LOG(xma)
  END SELECT
END FUNCTION QWGTS
