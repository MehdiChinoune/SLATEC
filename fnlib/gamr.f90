!DECK GAMR
FUNCTION GAMR(X)
  IMPLICIT NONE
  REAL alngx, GAMMA, GAMR, sgngx, X
  INTEGER irold
  !***BEGIN PROLOGUE  GAMR
  !***PURPOSE  Compute the reciprocal of the Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      SINGLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
  !***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! GAMR is a single precision function that evaluates the reciprocal
  ! of the gamma function for single precision argument X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ALGAMS, GAMMA, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  GAMR
  EXTERNAL GAMMA
  !***FIRST EXECUTABLE STATEMENT  GAMR
  GAMR = 0.0
  IF ( X<=0.0.AND.AINT(X)==X ) RETURN
  !
  CALL XGETF(irold)
  CALL XSETF(1)
  IF ( ABS(X)>10.0 ) THEN
    !
    CALL ALGAMS(X,alngx,sgngx)
    CALL XERCLR
    CALL XSETF(irold)
    GAMR = sgngx*EXP(-alngx)
    GOTO 99999
  ENDIF
  GAMR = 1.0/GAMMA(X)
  CALL XERCLR
  CALL XSETF(irold)
  RETURN
  !
  99999 CONTINUE
  END FUNCTION GAMR
