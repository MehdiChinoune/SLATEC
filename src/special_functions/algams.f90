!** ALGAMS
SUBROUTINE ALGAMS(X,Algam,Sgngam)
  !> Compute the logarithm of the absolute value of the Gamma
  !            function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      SINGLE PRECISION (ALGAMS-S, DLGAMS-D)
  !***
  ! **Keywords:**  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
  !             FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluates the logarithm of the absolute value of the gamma
  ! function.
  !     X           - input argument
  !     ALGAM       - result
  !     SGNGAM      - is set to the sign of GAMMA(X) and will
  !                   be returned at +1.0 or -1.0.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALNGAM

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  !* FIRST EXECUTABLE STATEMENT  ALGAMS
  REAL(SP) :: Algam, Sgngam, X
  INTEGER :: i
  Algam = LOG_GAMMA(X)
  Sgngam = 1.0
  IF( X>0.0 ) RETURN
  !
  i = INT( MOD(-AINT(X),2.0) + 0.1 )
  IF( i==0 ) Sgngam = -1.0
  !
END SUBROUTINE ALGAMS
