!** INITS
INTEGER PURE FUNCTION INITS(Os,Eta)
  !> Determine the number of terms needed in an orthogonal
  !  polynomial series so that it meets a specified accuracy.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C3A2
  !***
  ! **Type:**      SINGLE PRECISION (INITS-S, INITDS-D)
  !***
  ! **Keywords:**  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
  !             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891115  Modified error message.  (WRB)
  !   891115  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  REAL(SP), INTENT(IN) :: Eta
  REAL(SP), INTENT(IN) :: Os(:)
  INTEGER :: i, ii, nos
  REAL(SP) :: err
  !* FIRST EXECUTABLE STATEMENT  INITS
  nos = SIZE( Os )
  !
  err = 0._SP
  DO ii = 1, nos
    i = nos + 1 - ii
    err = err + ABS(Os(i))
    IF( err>Eta ) EXIT
  END DO
  !
  !IF( i==nos ) CALL XERMSG('INITS : Chebyshev series too short for specified accuracy',1,1)
  INITS = i
  !
END FUNCTION INITS
