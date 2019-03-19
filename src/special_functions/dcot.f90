!** DCOT
REAL(8) FUNCTION DCOT(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the cotangent.
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
  
  INTEGER ifn, INITDS, nterms
  REAL(8) :: X, cotcs(15), ainty, ainty2, pi2rec, sqeps, xmax, &
    xmin, xsml, y, yrem, prodbg, DCSEVL, D1MACH
  LOGICAL first
  SAVE cotcs, pi2rec, nterms, xmax, xsml, xmin, sqeps, first
  DATA cotcs(1)/ + .240259160982956302509553617744970D+0/
  DATA cotcs(2)/ - .165330316015002278454746025255758D-1/
  DATA cotcs(3)/ - .429983919317240189356476228239895D-4/
  DATA cotcs(4)/ - .159283223327541046023490851122445D-6/
  DATA cotcs(5)/ - .619109313512934872588620579343187D-9/
  DATA cotcs(6)/ - .243019741507264604331702590579575D-11/
  DATA cotcs(7)/ - .956093675880008098427062083100000D-14/
  DATA cotcs(8)/ - .376353798194580580416291539706666D-16/
  DATA cotcs(9)/ - .148166574646746578852176794666666D-18/
  DATA cotcs(10)/ - .583335658903666579477984000000000D-21/
  DATA cotcs(11)/ - .229662646964645773928533333333333D-23/
  DATA cotcs(12)/ - .904197057307483326719999999999999D-26/
  DATA cotcs(13)/ - .355988551920600064000000000000000D-28/
  DATA cotcs(14)/ - .140155139824298666666666666666666D-30/
  DATA cotcs(15)/ - .551800436872533333333333333333333D-33/
  DATA pi2rec/.011619772367581343075535053490057D0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  DCOT
  IF ( first ) THEN
    nterms = INITDS(cotcs,15,0.1*REAL(D1MACH(3)))
    xmax = 1.0D0/D1MACH(4)
    xsml = SQRT(3.0D0*D1MACH(3))
    xmin = EXP(MAX(LOG(D1MACH(1)),-LOG(D1MACH(2)))+0.01D0)
    sqeps = SQRT(D1MACH(4))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y<xmin ) CALL XERMSG('SLATEC','DCOT',&
    'ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS',2,2)
  IF ( y>xmax ) CALL XERMSG('SLATEC','DCOT',&
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
  IF ( ifn==1 ) y = 1.0D0 - y
  !
  IF ( ABS(X)>0.5D0.AND.y<ABS(X)*sqeps ) CALL XERMSG('SLATEC','DCOT',&
    'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI (N.NE.0)',1,1)
  !
  IF ( y<=0.25D0 ) THEN
    DCOT = 1.0D0/X
    IF ( y>xsml ) DCOT = (0.5D0+DCSEVL(32.0D0*y*y-1.D0,cotcs,nterms))/y
    !
  ELSEIF ( y>0.5D0 ) THEN
    !
    DCOT = (0.5D0+DCSEVL(2.D0*y*y-1.D0,cotcs,nterms))/(.25D0*y)
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
  ELSE
    DCOT = (0.5D0+DCSEVL(8.D0*y*y-1.D0,cotcs,nterms))/(0.5D0*y)
    DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
  ENDIF
  !
  IF ( X/=0.D0 ) DCOT = SIGN(DCOT,X)
  IF ( ifn==1 ) DCOT = -DCOT
  !
END FUNCTION DCOT
