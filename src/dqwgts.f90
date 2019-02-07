!*==DQWGTS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQWGTS
REAL(8) FUNCTION DQWGTS(X,A,B,Alfa,Beta,Integr)
  IMPLICIT NONE
  !*--DQWGTS5
  !***BEGIN PROLOGUE  DQWGTS
  !***SUBSIDIARY
  !***PURPOSE  This function subprogram is used together with the
  !            routine DQAWS and defines the WEIGHT function.
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (QWGTS-S, DQWGTS-D)
  !***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
  !             WEIGHT FUNCTION
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***SEE ALSO  DQK15W
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DQWGTS
  !
  REAL(8) :: A, Alfa, B, Beta, bmx, X, xma
  INTEGER Integr
  !***FIRST EXECUTABLE STATEMENT  DQWGTS
  xma = X - A
  bmx = B - X
  DQWGTS = xma**Alfa*bmx**Beta
  SELECT CASE (Integr)
    CASE (1)
    CASE (3)
      DQWGTS = DQWGTS*LOG(bmx)
    CASE (4)
      DQWGTS = DQWGTS*LOG(xma)*LOG(bmx)
    CASE DEFAULT
      DQWGTS = DQWGTS*LOG(xma)
  END SELECT
END FUNCTION DQWGTS
