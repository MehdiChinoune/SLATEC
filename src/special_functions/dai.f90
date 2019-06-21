!** DAI
REAL(DP) FUNCTION DAI(X)
  !> Evaluate the Airy function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      DOUBLE PRECISION (AI-S, DAI-D)
  !***
  ! **Keywords:**  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DAI(X) calculates the double precision Airy function for double
  ! precision argument X.
  !
  ! Series for AIF        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   8.37E-33
  !                                         log weighted error  32.08
  !                               significant figures required  30.87
  !                                    decimal places required  32.63
  !
  ! Series for AIG        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   7.47E-34
  !                                         log weighted error  33.13
  !                               significant figures required  31.50
  !                                    decimal places required  33.68
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9AIMP, DAIE, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, D1MACH
  REAL(DP) :: X
  REAL(DP) :: theta, xm, z
  INTEGER, SAVE :: naif, naig
  REAL(DP), PARAMETER :: x3sml = D1MACH(3)**0.3334_DP, &
    xmaxt = (-1.5_DP*LOG(D1MACH(1)))**0.6667_DP, &
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4._DP*SQRT(xmaxt)+1._DP) - 0.01_DP
  REAL(DP), PARAMETER :: aifcs(13) = [ -.37971358496669997496197089469414E-1_DP, &
    +.59191888537263638574319728013777E-1_DP, +.98629280577279975365603891044060E-3_DP, &
    +.68488438190765667554854830182412E-5_DP, +.25942025962194713019489279081403E-7_DP, &
    +.61766127740813750329445749697236E-10_DP, +.10092454172466117901429556224601E-12_DP, &
    +.12014792511179938141288033225333E-15_DP, +.10882945588716991878525295466666E-18_DP, &
    +.77513772196684887039238400000000E-22_DP, +.44548112037175638391466666666666E-25_DP, &
    +.21092845231692343466666666666666E-28_DP, +.83701735910741333333333333333333E-32_DP ]
  REAL(DP), PARAMETER :: aigcs(13) = [ +.18152365581161273011556209957864E-1_DP, &
    +.21572563166010755534030638819968E-1_DP, +.25678356987483249659052428090133E-3_DP, &
    +.14265214119792403898829496921721E-5_DP, +.45721149200180426070434097558191E-8_DP, &
    +.95251708435647098607392278840592E-11_DP, +.13925634605771399051150420686190E-13_DP, &
    +.15070999142762379592306991138666E-16_DP, +.12559148312567778822703205333333E-19_DP, &
    +.83063073770821340343829333333333E-23_DP, +.44657538493718567445333333333333E-26_DP, &
    +.19900855034518869333333333333333E-29_DP, +.74702885256533333333333333333333E-33_DP ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DAI
  IF( first ) THEN
    naif = INITDS(aifcs,13,0.1_SP*D1MACH(3))
    naig = INITDS(aigcs,13,0.1_SP*D1MACH(3))
    first = .FALSE.
  END IF
  !
  IF( X<(-1._DP) ) THEN
    CALL D9AIMP(X,xm,theta)
    DAI = xm*COS(theta)
    RETURN
    !
  ELSEIF( X<=1._DP ) THEN
    z = 0._DP
    IF( ABS(X)>x3sml ) z = X**3
    DAI = 0.375_DP + (DCSEVL(z,aifcs,naif)-X*(0.25_DP+DCSEVL(z,aigcs,naig)))
    RETURN
    !
  ELSEIF( X>xmax ) THEN
    !
    DAI = 0._DP
    CALL XERMSG('DAI','X SO BIG AI UNDERFLOWS',1,1)
    RETURN
  END IF
  DAI = DAIE(X)*EXP(-2._DP*X*SQRT(X)/3._DP)
  RETURN
END FUNCTION DAI
