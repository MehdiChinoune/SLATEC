!*==DAI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DAI
REAL(8) FUNCTION DAI(X)
  IMPLICIT NONE
  !*--DAI5
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS, naif, naig
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DAI
  !***PURPOSE  Evaluate the Airy function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      DOUBLE PRECISION (AI-S, DAI-D)
  !***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9AIMP, DAIE, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  !***END PROLOGUE  DAI
  REAL(8) :: X, aifcs(13), aigcs(13), theta, xm, xmax, x3sml, &
    z, D1MACH, DCSEVL, DAIE, xmaxt
  LOGICAL first
  SAVE aifcs, aigcs, naif, naig, x3sml, xmax, first
  DATA aifcs(1)/ - .37971358496669997496197089469414D-1/
  DATA aifcs(2)/ + .59191888537263638574319728013777D-1/
  DATA aifcs(3)/ + .98629280577279975365603891044060D-3/
  DATA aifcs(4)/ + .68488438190765667554854830182412D-5/
  DATA aifcs(5)/ + .25942025962194713019489279081403D-7/
  DATA aifcs(6)/ + .61766127740813750329445749697236D-10/
  DATA aifcs(7)/ + .10092454172466117901429556224601D-12/
  DATA aifcs(8)/ + .12014792511179938141288033225333D-15/
  DATA aifcs(9)/ + .10882945588716991878525295466666D-18/
  DATA aifcs(10)/ + .77513772196684887039238400000000D-22/
  DATA aifcs(11)/ + .44548112037175638391466666666666D-25/
  DATA aifcs(12)/ + .21092845231692343466666666666666D-28/
  DATA aifcs(13)/ + .83701735910741333333333333333333D-32/
  DATA aigcs(1)/ + .18152365581161273011556209957864D-1/
  DATA aigcs(2)/ + .21572563166010755534030638819968D-1/
  DATA aigcs(3)/ + .25678356987483249659052428090133D-3/
  DATA aigcs(4)/ + .14265214119792403898829496921721D-5/
  DATA aigcs(5)/ + .45721149200180426070434097558191D-8/
  DATA aigcs(6)/ + .95251708435647098607392278840592D-11/
  DATA aigcs(7)/ + .13925634605771399051150420686190D-13/
  DATA aigcs(8)/ + .15070999142762379592306991138666D-16/
  DATA aigcs(9)/ + .12559148312567778822703205333333D-19/
  DATA aigcs(10)/ + .83063073770821340343829333333333D-23/
  DATA aigcs(11)/ + .44657538493718567445333333333333D-26/
  DATA aigcs(12)/ + .19900855034518869333333333333333D-29/
  DATA aigcs(13)/ + .74702885256533333333333333333333D-33/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DAI
  IF ( first ) THEN
    naif = INITDS(aifcs,13,0.1*REAL(D1MACH(3)))
    naig = INITDS(aigcs,13,0.1*REAL(D1MACH(3)))
    !
    x3sml = D1MACH(3)**0.3334D0
    xmaxt = (-1.5D0*LOG(D1MACH(1)))**0.6667D0
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4.0D0*SQRT(xmaxt)+1.0D0) - 0.01D0
  ENDIF
  first = .FALSE.
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
    GOTO 99999
  ENDIF
  DAI = DAIE(X)*EXP(-2.0D0*X*SQRT(X)/3.0D0)
  RETURN
  !
  99999 CONTINUE
  END FUNCTION DAI
