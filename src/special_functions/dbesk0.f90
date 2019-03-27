!** DBESK0
REAL(8) FUNCTION DBESK0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the modified (hyperbolic) Bessel function of the
  !            third kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      DOUBLE PRECISION (BESK0-S, DBESK0-D)
  !***
  ! **Keywords:**  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESK0(X) calculates the double precision modified (hyperbolic)
  ! Bessel function of the third kind of order zero for double
  ! precision argument X.  The argument must be greater than zero
  ! but not so large that the result underflows.
  !
  ! Series for BK0        on the interval  0.          to  4.00000E+00
  !                                        with weighted error   3.08E-33
  !                                         log weighted error  32.51
  !                               significant figures required  32.05
  !                                    decimal places required  33.11
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  INTEGER INITDS, ntk0
  REAL(8) :: X, xmax, xmaxt, xsml, y, D1MACH, DCSEVL, DBESI0, DBSK0E
  SAVE ntk0, xsml, xmax
  REAL(8), PARAMETER :: bk0cs(16) = [ -.353273932339027687201140060063153D-1, &
    +.344289899924628486886344927529213D+0, +.359799365153615016265721303687231D-1, &
    +.126461541144692592338479508673447D-2, +.228621210311945178608269830297585D-4, &
    +.253479107902614945730790013428354D-6, +.190451637722020885897214059381366D-8, &
    +.103496952576336245851008317853089D-10, +.425981614279108257652445327170133D-13, &
    +.137446543588075089694238325440000D-15, +.357089652850837359099688597333333D-18, &
    +.763164366011643737667498666666666D-21, +.136542498844078185908053333333333D-23, &
    +.207527526690666808319999999999999D-26, +.271281421807298560000000000000000D-29, &
    +.308259388791466666666666666666666D-32 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBESK0
  IF ( first ) THEN
    ntk0 = INITDS(bk0cs,16,0.1*REAL(D1MACH(3)))
    xsml = SQRT(4.0D0*D1MACH(3))
    xmaxt = -LOG(D1MACH(1))
    xmax = xmaxt - 0.5D0*xmaxt*LOG(xmaxt)/(xmaxt+0.5D0)
    first = .FALSE.
  ENDIF
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBESK0','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2.0D0 ) THEN
    !
    DBESK0 = 0.D0
    IF ( X>xmax ) CALL XERMSG('SLATEC','DBESK0','X SO BIG K0 UNDERFLOWS',1,1)
    IF ( X>xmax ) RETURN
    !
    DBESK0 = EXP(-X)*DBSK0E(X)
    RETURN
  ENDIF
  !
  y = 0.D0
  IF ( X>xsml ) y = X*X
  DBESK0 = -LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + DCSEVL(.5D0*y-1.D0,bk0cs,ntk0)
  RETURN
END FUNCTION DBESK0
