!** DCOT
REAL(DP) FUNCTION DCOT(X)
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
  USE service, ONLY : XERMSG, D1MACH
  REAL(DP) :: X
  INTEGER :: ifn
  REAL(DP) :: ainty, ainty2, y, yrem, prodbg
  INTEGER, SAVE :: nterms
  REAL(DP), PARAMETER :: xmax = 1.0D0/D1MACH(4), xsml = SQRT(3.0D0*D1MACH(3)), &
    xmin = EXP(MAX(LOG(D1MACH(1)),-LOG(D1MACH(2)))+0.01D0), sqeps = SQRT(D1MACH(4))
  REAL(DP), PARAMETER :: cotcs(15) = [ +.240259160982956302509553617744970D+0, &
    -.165330316015002278454746025255758D-1, -.429983919317240189356476228239895D-4, &
    -.159283223327541046023490851122445D-6, -.619109313512934872588620579343187D-9, &
    -.243019741507264604331702590579575D-11, -.956093675880008098427062083100000D-14, &
    -.376353798194580580416291539706666D-16, -.148166574646746578852176794666666D-18, &
    -.583335658903666579477984000000000D-21, -.229662646964645773928533333333333D-23, &
    -.904197057307483326719999999999999D-26, -.355988551920600064000000000000000D-28, &
    -.140155139824298666666666666666666D-30, -.551800436872533333333333333333333D-33 ]
  REAL(DP), PARAMETER :: pi2rec = .011619772367581343075535053490057D0
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DCOT
  IF( first ) THEN
    nterms = INITDS(cotcs,15,0.1*D1MACH(3))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF( y<xmin ) CALL XERMSG('DCOT',&
    'ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS',2,2)
  IF( y>xmax ) CALL XERMSG('DCOT',&
    'NO PRECISION BECAUSE ABS(X) IS TOO BIG',3,2)
  !
  ! CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
  ! = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
  ! = AINT(.625*Y) + AINT(Z) + REM(Z)
  !
  ainty = AINT(y)
  yrem = y - ainty
  prodbg = 0.625D0*ainty
  ainty = AINT(prodbg)
  y = (prodbg-ainty) + 0.625D0*yrem + pi2rec*y
  ainty2 = AINT(y)
  ainty = ainty + ainty2
  y = y - ainty2
  !
  ifn = INT( MOD(ainty,2.0D0) )
  IF( ifn==1 ) y = 1.0D0 - y
  !
  IF( ABS(X)>0.5D0 .AND. y<ABS(X)*sqeps ) CALL XERMSG('DCOT',&
    'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI (N/=0)',1,1)
  !
  IF( y<=0.25D0 ) THEN
    DCOT = 1.0D0/X
    IF( y>xsml ) DCOT = (0.5D0+DCSEVL(32.0D0*y*y-1.D0,cotcs,nterms))/y
    !
  ELSEIF( y>0.5D0 ) THEN
    !
    DCOT = (0.5D0+DCSEVL(2.D0*y*y-1.D0,cotcs,nterms))/(.25D0*y)
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
  ELSE
    DCOT = (0.5D0+DCSEVL(8.D0*y*y-1.D0,cotcs,nterms))/(0.5D0*y)
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
  END IF
  !
  IF( X/=0.D0 ) DCOT = SIGN(DCOT,X)
  IF( ifn==1 ) DCOT = -DCOT
  !
END FUNCTION DCOT
