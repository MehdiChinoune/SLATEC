!** GAMLIM
SUBROUTINE GAMLIM(Xmin,Xmax)
  IMPLICIT NONE
  !>
  !***
  !  Compute the minimum and maximum bounds for the argument in
  !            the Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A, R2
  !***
  ! **Type:**      SINGLE PRECISION (GAMLIM-S, DGAMLM-D)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Calculate the minimum and maximum legal bounds for X in GAMMA(X).
  ! XMIN and XMAX are not the only bounds, but they are the only non-
  ! trivial ones to calculate.
  !
  !             Output Arguments --
  ! XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of
  !        X might result in underflow.
  ! XMAX   maximum legal value of X in GAMMA(X).  Any larger value will
  !        cause overflow.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  
  REAL alnbig, alnsml, R1MACH, xln, Xmax, Xmin, xold
  INTEGER i
  !* FIRST EXECUTABLE STATEMENT  GAMLIM
  alnsml = LOG(R1MACH(1))
  Xmin = -alnsml
  DO i = 1, 10
    xold = Xmin
    xln = LOG(Xmin)
    Xmin = Xmin - Xmin*((Xmin+0.5)*xln-Xmin-0.2258+alnsml)/(Xmin*xln+0.5)
    IF ( ABS(Xmin-xold)<0.005 ) GOTO 100
  END DO
  CALL XERMSG('SLATEC','GAMLIM','UNABLE TO FIND XMIN',1,2)
  !
  100  Xmin = -Xmin + 0.01
  !
  alnbig = LOG(R1MACH(2))
  Xmax = alnbig
  DO i = 1, 10
    xold = Xmax
    xln = LOG(Xmax)
    Xmax = Xmax - Xmax*((Xmax-0.5)*xln-Xmax+0.9189-alnbig)/(Xmax*xln-0.5)
    IF ( ABS(Xmax-xold)<0.005 ) GOTO 200
  END DO
  CALL XERMSG('SLATEC','GAMLIM','UNABLE TO FIND XMAX',2,2)
  !
  200  Xmax = Xmax - 0.01
  Xmin = MAX(Xmin,-Xmax+1.)
  !
END SUBROUTINE GAMLIM
