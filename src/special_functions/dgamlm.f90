!** DGAMLM
SUBROUTINE DGAMLM(Xmin,Xmax)
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
  ! **Type:**      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Calculate the minimum and maximum legal bounds for X in gamma(X).
  ! XMIN and XMAX are not the only bounds, but they are the only non-
  ! trivial ones to calculate.
  !
  !             Output Arguments --
  ! XMIN   double precision minimum legal value of X in gamma(X).  Any
  !        smaller value of X might result in underflow.
  ! XMAX   double precision maximum legal value of X in gamma(X).  Any
  !        larger value of X might cause overflow.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  
  INTEGER i
  REAL(8) :: Xmin, Xmax, alnbig, alnsml, xln, xold, D1MACH
  !* FIRST EXECUTABLE STATEMENT  DGAMLM
  alnsml = LOG(D1MACH(1))
  Xmin = -alnsml
  DO i = 1, 10
    xold = Xmin
    xln = LOG(Xmin)
    Xmin = Xmin - Xmin*((Xmin+0.5D0)*xln-Xmin-0.2258D0+alnsml)/(Xmin*xln+0.5D0)
    IF ( ABS(Xmin-xold)<0.005D0 ) GOTO 100
  ENDDO
  CALL XERMSG('SLATEC','DGAMLM','UNABLE TO FIND XMIN',1,2)
  !
  100  Xmin = -Xmin + 0.01D0
  !
  alnbig = LOG(D1MACH(2))
  Xmax = alnbig
  DO i = 1, 10
    xold = Xmax
    xln = LOG(Xmax)
    Xmax = Xmax - Xmax*((Xmax-0.5D0)*xln-Xmax+0.9189D0-alnbig)/(Xmax*xln-0.5D0)
    IF ( ABS(Xmax-xold)<0.005D0 ) GOTO 200
  ENDDO
  CALL XERMSG('SLATEC','DGAMLM','UNABLE TO FIND XMAX',2,2)
  !
  200  Xmax = Xmax - 0.01D0
  Xmin = MAX(Xmin,-Xmax+1.D0)
  !
END SUBROUTINE DGAMLM
