!*==INITDS.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK INITDS
FUNCTION INITDS(Os,Nos,Eta)
  IMPLICIT NONE
  !*--INITDS5
  !*** Start of declarations inserted by SPAG
  REAL err , Eta
  INTEGER i , ii , INITDS , Nos
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  INITDS
  !***PURPOSE  Determine the number of terms needed in an orthogonal
  !            polynomial series so that it meets a specified accuracy.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C3A2
  !***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
  !***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
  !             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !  Initialize the orthogonal series, represented by the array OS, so
  !  that INITDS is the number of terms needed to insure the error is no
  !  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
  !  machine precision.
  !
  !             Input Arguments --
  !   OS     double precision array of NOS coefficients in an orthogonal
  !          series.
  !   NOS    number of coefficients in OS.
  !   ETA    single precision scalar containing requested accuracy of
  !          series.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891115  Modified error message.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  INITDS
  DOUBLE PRECISION Os(*)
  !***FIRST EXECUTABLE STATEMENT  INITDS
  IF ( Nos<1 ) CALL XERMSG('SLATEC','INITDS',&
    'Number of coefficients is less than 1',2,1)
  !
  err = 0.
  DO ii = 1 , Nos
    i = Nos + 1 - ii
    err = err + ABS(REAL(Os(i)))
    IF ( err>Eta ) EXIT
  ENDDO
  !
  IF ( i==Nos ) CALL XERMSG('SLATEC','INITDS',&
    'Chebyshev series too short for specified accuracy'&
    ,1,1)
  INITDS = i
  !
END FUNCTION INITDS
