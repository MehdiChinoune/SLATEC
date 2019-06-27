!** COT
REAL(SP) ELEMENTAL FUNCTION COT(X)
  !> Compute the cotangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      SINGLE PRECISION (COT-S, DCOT-D, CCOT-C)
  !***
  ! **Keywords:**  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! COT(X) calculates the cotangent of the real argument X.  X is in
  ! units of radians.
  !
  ! Series for COT        on the interval  0.          to  6.25000D-02
  !                                        with weighted error   3.76E-17
  !                                         log weighted error  16.42
  !                               significant figures required  15.51
  !                                    decimal places required  16.88
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  INTEGER :: ifn
  REAL(SP) :: ainty, ainty2, prodbg, y, yrem
  INTEGER, PARAMETER :: nterms = 4
  REAL(SP), PARAMETER :: xmax = 1._SP/R1MACH(4), xsml = SQRT(3._SP*R1MACH(3)), &
    xmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+0.01_SP), sqeps = SQRT(R1MACH(4))
  REAL(SP), PARAMETER :: cotcs(8) = [ .24025916098295630_SP,-.016533031601500228_SP, &
    -.000042998391931724_SP,-.000000159283223327_SP,-.000000000619109313_SP, &
    -.000000000002430197_SP,-.000000000000009560_SP,-.000000000000000037_SP ]
  REAL(SP), PARAMETER :: pi2rec = .0116197723675813430_SP
  !* FIRST EXECUTABLE STATEMENT  COT
  ! nterms = INITS(cotcs,0.1_SP*R1MACH(3))
  !
  y = ABS(X)
  IF( y<xmin ) ERROR STOP 'COT : ABS(X) IS ZERO OR SO SMALL COT OVERFLOWS'
  IF( y>xmax ) ERROR STOP 'COT : NO PRECISION BECAUSE ABS(X) IS TOO BIG'
  !
  ! CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
  ! = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
  ! = AINT(.625*Y) + AINT(Z) + REM(Z)
  !
  ainty = AINT(y)
  yrem = y - ainty
  prodbg = 0.625_SP*ainty
  ainty = AINT(prodbg)
  y = (prodbg-ainty) + 0.625_SP*yrem + y*pi2rec
  ainty2 = AINT(y)
  ainty = ainty + ainty2
  y = y - ainty2
  !
  ifn = INT( MOD(ainty,2._SP) )
  IF( ifn==1 ) y = 1._SP - y
  !
  ! IF( ABS(X)>0.5_SP .AND. y<ABS(X)*sqeps ) CALL XERMSG('COT',&
    ! 'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI (N/=0)',1,1)
  !
  IF( y<=0.25 ) THEN
    COT = 1._SP/X
    IF( y>xsml ) COT = (0.5_SP+CSEVL(32._SP*y*y-1._SP,cotcs(1:nterms)))/y
    !
  ELSEIF( y>0.5_SP ) THEN
    !
    COT = (0.5_SP+CSEVL(2._SP*y*y-1._SP,cotcs(1:nterms)))/(0.25_SP*y)
    COT = (COT**2-1._SP)*0.5_SP/COT
    COT = (COT**2-1._SP)*0.5_SP/COT
  ELSE
    COT = (0.5_SP+CSEVL(8._SP*y*y-1._SP,cotcs(1:nterms)))/(0.5_SP*y)
    COT = (COT**2-1._SP)*0.5_SP/COT
  END IF
  !
  COT = SIGN(COT,X)
  IF( ifn==1 ) COT = -COT
  !
END FUNCTION COT