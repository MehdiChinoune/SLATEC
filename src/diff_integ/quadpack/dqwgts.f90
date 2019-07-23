!** DQWGTS
REAL(DP) PURE FUNCTION DQWGTS(X,A,B,Alfa,Beta,Integr)
  !> This function subprogram is used together with the routine
  !  DQAWS and defines the WEIGHT function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (QWGTS-S, DQWGTS-D)
  !***
  ! **Keywords:**  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
  !             WEIGHT FUNCTION
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **See also:**  DQK15W
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  INTEGER, INTENT(IN) :: Integr
  REAL(DP), INTENT(IN) :: A, Alfa, B, Beta, X
  !
  REAL(DP) :: bmx, xma
  !* FIRST EXECUTABLE STATEMENT  DQWGTS
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
  !
END FUNCTION DQWGTS