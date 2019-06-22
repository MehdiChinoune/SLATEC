!** DGAMR
REAL(DP) FUNCTION DGAMR(X)
  !> Compute the reciprocal of the Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      DOUBLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
  !***
  ! **Keywords:**  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DGAMR(X) calculates the double precision reciprocal of the
  ! complete Gamma function for double precision argument X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DLGAMS, XERCLR, XGETF, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : control_xer, num_xer
  INTEGER :: irold
  REAL(DP) :: X, alngx, sgngx
  !* FIRST EXECUTABLE STATEMENT  DGAMR
  DGAMR = 0._DP
  IF( X<=0._DP .AND. AINT(X)==X ) RETURN
  !
  irold = control_xer
  control_xer = 1
  IF( ABS(X)>10._DP ) THEN
    !
    CALL DLGAMS(X,alngx,sgngx)
    num_xer = 0
    control_xer = irold
    DGAMR = sgngx*EXP(-alngx)
    RETURN
  END IF
  DGAMR = 1._DP/GAMMA(X)
  num_xer = 0
  control_xer = irold
  RETURN
END FUNCTION DGAMR
