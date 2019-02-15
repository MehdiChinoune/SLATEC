!DECK DSOSFN
REAL(8) FUNCTION DSOSFN(X,K)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DSOSFN
  !***PURPOSE  Function evaluator for DSOS quick check.
  !***LIBRARY   SLATEC
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !     FUNCTION WHICH EVALUATES THE FUNCTIONS, ONE AT A TIME,
  !     FOR TEST PROGRAM USED IN QUICK CHECK OF DSOS.
  !
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DSOSFN
  INTEGER K
  REAL(8) :: X(2)
  !***FIRST EXECUTABLE STATEMENT  DSOSFN
  IF ( K==1 ) DSOSFN = 1.0D0 - X(1)
  IF ( K==2 ) DSOSFN = 1.0D1*(X(2)-X(1)**2)
END FUNCTION DSOSFN
