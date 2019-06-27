!** DCOT
REAL(DP) ELEMENTAL FUNCTION DCOT(X)
  !> Compute the cotangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      DOUBLE PRECISION (COT-S, DCOT-D, CCOT-C)
  !***
  ! **Keywords:**  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DCOT(X) calculates the double precision trigonometric cotangent
  ! for double precision argument X.  X is in units of radians.
  !
  ! Series for COT        on the interval  0.          to  6.25000E-02
  !                                        with weighted error   5.52E-34
  !                                         log weighted error  33.26
  !                               significant figures required  32.34
  !                                    decimal places required  33.85
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: X
  INTEGER :: ifn
  REAL(DP) :: ainty, ainty2, y, yrem, prodbg
  INTEGER, PARAMETER :: nterms = 8
  REAL(DP), PARAMETER :: xmax = 1._DP/D1MACH(4), xsml = SQRT(3._DP*D1MACH(3)), &
    xmin = EXP(MAX(LOG(D1MACH(1)),-LOG(D1MACH(2)))+0.01_DP), sqeps = SQRT(D1MACH(4))
  REAL(DP), PARAMETER :: cotcs(15) = [ +.240259160982956302509553617744970E+0_DP, &
    -.165330316015002278454746025255758E-1_DP, -.429983919317240189356476228239895E-4_DP, &
    -.159283223327541046023490851122445E-6_DP, -.619109313512934872588620579343187E-9_DP, &
    -.243019741507264604331702590579575E-11_DP, -.956093675880008098427062083100000E-14_DP, &
    -.376353798194580580416291539706666E-16_DP, -.148166574646746578852176794666666E-18_DP, &
    -.583335658903666579477984000000000E-21_DP, -.229662646964645773928533333333333E-23_DP, &
    -.904197057307483326719999999999999E-26_DP, -.355988551920600064000000000000000E-28_DP, &
    -.140155139824298666666666666666666E-30_DP, -.551800436872533333333333333333333E-33_DP ]
  REAL(DP), PARAMETER :: pi2rec = .011619772367581343075535053490057_DP
  !* FIRST EXECUTABLE STATEMENT  DCOT
  ! nterms = INITDS(cotcs,0.1_SP*D1MACH(3))
  !
  y = ABS(X)
  IF( y<xmin ) ERROR STOP 'DCOT : ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS'
  IF( y>xmax ) ERROR STOP 'DCOT : NO PRECISION BECAUSE ABS(X) IS TOO BIG'
  !
  ! CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
  ! = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
  ! = AINT(.625*Y) + AINT(Z) + REM(Z)
  !
  ainty = AINT(y)
  yrem = y - ainty
  prodbg = 0.625_DP*ainty
  ainty = AINT(prodbg)
  y = (prodbg-ainty) + 0.625_DP*yrem + pi2rec*y
  ainty2 = AINT(y)
  ainty = ainty + ainty2
  y = y - ainty2
  !
  ifn = INT( MOD(ainty,2._DP) )
  IF( ifn==1 ) y = 1._DP - y
  !
  ! IF( ABS(X)>0.5_DP .AND. y<ABS(X)*sqeps ) CALL XERMSG('DCOT',&
    ! 'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI (N/=0)',1,1)
  !
  IF( y<=xsml ) THEN
    DCOT = 1._DP/X
  ELSEIF( y<=0.25_DP ) THEN
    DCOT = (0.5_DP+DCSEVL(32._DP*y*y-1._DP,cotcs(1:nterms)))/y
  ELSEIF( y<=0.5_DP ) THEN
    DCOT = (0.5_DP+DCSEVL(8._DP*y*y-1._DP,cotcs(1:nterms)))/(0.5_DP*y)
    DCOT = (DCOT*DCOT-1._DP)*0.5_DP/DCOT
  ELSE
    DCOT = (0.5_DP+DCSEVL(2._DP*y*y-1._DP,cotcs(1:nterms)))/(.25_DP*y)
    DCOT = (DCOT*DCOT-1._DP)*0.5_DP/DCOT
    DCOT = (DCOT*DCOT-1._DP)*0.5_DP/DCOT
  END IF
  !
  DCOT = SIGN(DCOT,X)
  IF( ifn==1 ) DCOT = -DCOT
  !
END FUNCTION DCOT