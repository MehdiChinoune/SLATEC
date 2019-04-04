!** DAI
REAL(8) FUNCTION DAI(X)
  IMPLICIT NONE
  !>
  !***
  !  Evaluate the Airy function.
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

  INTEGER INITDS
  REAL(8) :: X, theta, xm, z, D1MACH, DCSEVL, DAIE, xmaxt
  INTEGER, SAVE :: naif, naig
  REAL(8), SAVE :: x3sml, xmax
  REAL(8), PARAMETER :: aifcs(13) = [ -.37971358496669997496197089469414D-1, &
    +.59191888537263638574319728013777D-1, +.98629280577279975365603891044060D-3, &
    +.68488438190765667554854830182412D-5, +.25942025962194713019489279081403D-7, &
    +.61766127740813750329445749697236D-10, +.10092454172466117901429556224601D-12, &
    +.12014792511179938141288033225333D-15, +.10882945588716991878525295466666D-18, &
    +.77513772196684887039238400000000D-22, +.44548112037175638391466666666666D-25, &
    +.21092845231692343466666666666666D-28, +.83701735910741333333333333333333D-32 ]
  REAL(8), PARAMETER :: aigcs(13) = [ +.18152365581161273011556209957864D-1, &
    +.21572563166010755534030638819968D-1, +.25678356987483249659052428090133D-3, &
    +.14265214119792403898829496921721D-5, +.45721149200180426070434097558191D-8, &
    +.95251708435647098607392278840592D-11, +.13925634605771399051150420686190D-13, &
    +.15070999142762379592306991138666D-16, +.12559148312567778822703205333333D-19, &
    +.83063073770821340343829333333333D-23, +.44657538493718567445333333333333D-26, &
    +.19900855034518869333333333333333D-29, +.74702885256533333333333333333333D-33 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DAI
  IF ( first ) THEN
    naif = INITDS(aifcs,13,0.1*REAL(D1MACH(3)))
    naig = INITDS(aigcs,13,0.1*REAL(D1MACH(3)))
    !
    x3sml = D1MACH(3)**0.3334D0
    xmaxt = (-1.5D0*LOG(D1MACH(1)))**0.6667D0
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4.0D0*SQRT(xmaxt)+1.0D0) - 0.01D0
    first = .FALSE.
  END IF
  !
  IF ( X<(-1.D0) ) THEN
    CALL D9AIMP(X,xm,theta)
    DAI = xm*COS(theta)
    RETURN
    !
  ELSEIF ( X<=1.0D0 ) THEN
    z = 0.0D0
    IF ( ABS(X)>x3sml ) z = X**3
    DAI = 0.375D0 + (DCSEVL(z,aifcs,naif)-X*(0.25D0+DCSEVL(z,aigcs,naig)))
    RETURN
    !
  ELSEIF ( X>xmax ) THEN
    !
    DAI = 0.0D0
    CALL XERMSG('SLATEC','DAI','X SO BIG AI UNDERFLOWS',1,1)
    RETURN
  END IF
  DAI = DAIE(X)*EXP(-2.0D0*X*SQRT(X)/3.0D0)
  RETURN
END FUNCTION DAI
