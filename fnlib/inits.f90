!DECK INITS
FUNCTION INITS(Os,Nos,Eta)
  IMPLICIT NONE
  REAL err, Eta
  INTEGER i, ii, INITS, Nos
  !***BEGIN PROLOGUE  INITS
  !***PURPOSE  Determine the number of terms needed in an orthogonal
  !            polynomial series so that it meets a specified accuracy.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C3A2
  !***TYPE      SINGLE PRECISION (INITS-S, INITDS-D)
  !***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
  !             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !  Initialize the orthogonal series, represented by the array OS, so
  !  that INITS is the number of terms needed to insure the error is no
  !  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
  !  machine precision.
  !
  !             Input Arguments --
  !   OS     single precision array of NOS coefficients in an orthogonal
  !          series.
  !   NOS    number of coefficients in OS.
  !   ETA    single precision scalar containing requested accuracy of
  !          series.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891115  Modified error message.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  INITS
  REAL Os(*)
  !***FIRST EXECUTABLE STATEMENT  INITS
  IF ( Nos<1 ) CALL XERMSG('SLATEC','INITS',&
    'Number of coefficients is less than 1',2,1)
  !
  err = 0.
  DO ii = 1, Nos
    i = Nos + 1 - ii
    err = err + ABS(Os(i))
    IF ( err>Eta ) EXIT
  ENDDO
  !
  IF ( i==Nos ) CALL XERMSG('SLATEC','INITS',&
    'Chebyshev series too short for specified accuracy'&
    ,1,1)
  INITS = i
  !
END FUNCTION INITS
