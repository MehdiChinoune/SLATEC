!** QWGTF
REAL(SP) PURE FUNCTION QWGTF(X,Omega,P2,P3,P4,Integr)
  !> This function subprogram is used together with the routine
  !  QAWF and defines the WEIGHT function.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (QWGTF-S, DQWGTF-D)
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
  ! **See also:**  QK15W
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   810101  DATE WRITTEN
  !   830518  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  INTEGER, INTENT(IN) :: Integr
  REAL(SP), INTENT(IN) :: Omega, P2, P3, P4, X
  !
  REAL(SP) :: omx
  !* FIRST EXECUTABLE STATEMENT  QWGTF
  omx = Omega*X
  IF( Integr==2 ) THEN
    QWGTF = SIN(omx)
  ELSE
    QWGTF = COS(omx)
  END IF
  !
END FUNCTION QWGTF