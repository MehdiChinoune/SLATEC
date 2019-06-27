!** DBESK0
REAL(DP) ELEMENTAL FUNCTION DBESK0(X)
  !> Compute the modified (hyperbolic) Bessel function of the third kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      DOUBLE PRECISION (BESK0-S, DBESK0-D)
  !***
  ! **Keywords:**  FNLIB, HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
  !                ORDER ZERO, SPECIAL FUNCTIONS, THIRD KIND
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
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: X
  REAL(DP) :: y
  INTEGER, PARAMETER :: ntk0 = 10
  REAL(DP), PARAMETER :: xsml = SQRT(4._DP*D1MACH(3)), xmaxt = -LOG(D1MACH(1)), &
    xmax = xmaxt - 0.5_DP*xmaxt*LOG(xmaxt)/(xmaxt+0.5_DP)
  REAL(DP), PARAMETER :: bk0cs(16) = [ -.353273932339027687201140060063153E-1_DP, &
    +.344289899924628486886344927529213E+0_DP, +.359799365153615016265721303687231E-1_DP, &
    +.126461541144692592338479508673447E-2_DP, +.228621210311945178608269830297585E-4_DP, &
    +.253479107902614945730790013428354E-6_DP, +.190451637722020885897214059381366E-8_DP, &
    +.103496952576336245851008317853089E-10_DP, +.425981614279108257652445327170133E-13_DP, &
    +.137446543588075089694238325440000E-15_DP, +.357089652850837359099688597333333E-18_DP, &
    +.763164366011643737667498666666666E-21_DP, +.136542498844078185908053333333333E-23_DP, &
    +.207527526690666808319999999999999E-26_DP, +.271281421807298560000000000000000E-29_DP, &
    +.308259388791466666666666666666666E-32_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBESK0
  ! ntk0 = INITDS(bk0cs,0.1_SP*D1MACH(3))
  !
  IF( X<=0._DP ) THEN
    ERROR STOP 'DBESK0 : X IS ZERO OR NEGATIVE'
  ELSEIF( X<=2._DP ) THEN
    y = 0._DP
    IF( X>xsml ) y = X*X
    DBESK0 = -LOG(0.5_DP*X)*DBESI0(X) - 0.25_DP + DCSEVL(.5_DP*y-1._DP,bk0cs(1:ntk0))
  ELSEIF( X<=xmax ) THEN
    DBESK0 = EXP(-X)*DBSK0E(X)
  ELSE
    DBESK0 = 0._DP
    ! IF( X>xmax ) CALL XERMSG('DBESK0','X SO BIG K0 UNDERFLOWS',1,1)
  END IF

  RETURN
END FUNCTION DBESK0