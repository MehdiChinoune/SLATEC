!*==DGAMR.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DGAMR
REAL(8) FUNCTION DGAMR(X)
  IMPLICIT NONE
  !*--DGAMR5
  !*** Start of declarations inserted by SPAG
  INTEGER irold
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DGAMR
  !***PURPOSE  Compute the reciprocal of the Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      DOUBLE PRECISION (GAMR-S, DGAMR-D, CGAMR-C)
  !***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DGAMR(X) calculates the double precision reciprocal of the
  ! complete Gamma function for double precision argument X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DGAMMA, DLGAMS, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  DGAMR
  REAL(8) :: X , alngx , sgngx , DGAMMA
  EXTERNAL DGAMMA
  !***FIRST EXECUTABLE STATEMENT  DGAMR
  DGAMR = 0.0D0
  IF ( X<=0.0D0.AND.AINT(X)==X ) RETURN
  !
  CALL XGETF(irold)
  CALL XSETF(1)
  IF ( ABS(X)>10.0D0 ) THEN
    !
    CALL DLGAMS(X,alngx,sgngx)
    CALL XERCLR
    CALL XSETF(irold)
    DGAMR = sgngx*EXP(-alngx)
    GOTO 99999
  ENDIF
  DGAMR = 1.0D0/DGAMMA(X)
  CALL XERCLR
  CALL XSETF(irold)
  RETURN
  !
  99999 END FUNCTION DGAMR
