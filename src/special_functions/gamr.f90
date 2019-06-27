!** GAMR
REAL(SP) ELEMENTAL FUNCTION GAMR(X)
  !> Compute the reciprocal of the Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      SINGLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
  !***
  ! **Keywords:**  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! GAMR is a single precision function that evaluates the reciprocal
  ! of the gamma function for single precision argument X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALGAMS, XERCLR, XGETF, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: alngx, sgngx
  !* FIRST EXECUTABLE STATEMENT  GAMR

  IF( X<=0._SP .AND. AINT(X)==X ) THEN
    GAMR = 0._SP
  ELSEIF( ABS(X)>10._SP ) THEN
    CALL ALGAMS(X,alngx,sgngx)
    GAMR = sgngx*EXP(-alngx)
  ELSE
    GAMR = 1._SP/GAMMA(X)
  END IF

  RETURN
END FUNCTION GAMR