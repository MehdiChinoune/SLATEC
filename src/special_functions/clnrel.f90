!DECK CLNREL
COMPLEX FUNCTION CLNREL(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CLNREL
  !***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4B
  !***TYPE      COMPLEX (ALNREL-S, DLNREL-D, CLNREL-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CLNREL(Z) = LOG(1+Z) with relative error accuracy near Z = 0.
  ! Let   RHO = ABS(Z)  and
  !       R**2 = ABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 .
  ! Now if RHO is small we may evaluate CLNREL(Z) accurately by
  !       LOG(1+Z) = CMPLX  (LOG(R), CARG(1+Z))
  !                 = CMPLX  (0.5*LOG(R**2), CARG(1+Z))
  !                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z))
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ALNREL, CARG, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  CLNREL
  REAL ALNREL, CARG, R1MACH, rho, sqeps, x
  COMPLEX Z
  SAVE sqeps
  DATA sqeps/0.0/
  !***FIRST EXECUTABLE STATEMENT  CLNREL
  IF ( sqeps==0. ) sqeps = SQRT(R1MACH(4))
  !
  IF ( ABS(1.+Z)<sqeps ) CALL XERMSG('SLATEC','CLNREL',&
    'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR -1',1,1)
  !
  rho = ABS(Z)
  IF ( rho>0.375 ) CLNREL = LOG(1.0+Z)
  IF ( rho>0.375 ) RETURN
  !
  x = REAL(Z)
  CLNREL = CMPLX(0.5*ALNREL(2.*x+rho**2),CARG(1.0+Z))
  !
END FUNCTION CLNREL
