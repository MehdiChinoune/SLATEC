!** DQWGTF
REAL(DP) FUNCTION DQWGTF(X,Omega,P2,P3,P4,Integr)
  !>
  !  This function subprogram is used together with the
  !            routine DQAWF and defines the WEIGHT function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (QWGTF-S, DQWGTF-D)
  !***
  ! **Keywords:**  COS OR SIN IN WEIGHT FUNCTION
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
  REAL(DP) :: Omega, omx, P2, P3, P4, X
  INTEGER Integr
  !* FIRST EXECUTABLE STATEMENT  DQWGTF
  omx = Omega*X
  IF ( Integr==2 ) THEN
    DQWGTF = SIN(omx)
  ELSE
    DQWGTF = COS(omx)
  END IF
END FUNCTION DQWGTF
