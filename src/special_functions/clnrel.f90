!** CLNREL
COMPLEX(SP) ELEMENTAL FUNCTION CLNREL(Z)
  !> Evaluate ln(1+X) accurate in the sense of relative error.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      COMPLEX (ALNREL-S, DLNREL-D, CLNREL-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CLNREL(Z) = LOG(1+Z) with relative error accuracy near Z = 0.
  ! Let   RHO = ABS(Z)  and
  !       R**2 = ABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 .
  ! Now if RHO is small we may evaluate CLNREL(Z) accurately by
  !       LOG(1+Z) = CMPLX  (LOG(R), CARG(1+Z))
  !                 = CMPLX  (0.5*LOG(R**2), CARG(1+Z))
  !                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z))
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALNREL, CARG, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !
  COMPLEX(SP), INTENT(IN) :: Z
  !
  REAL(SP) :: rho, x
  !* FIRST EXECUTABLE STATEMENT  CLNREL
  !
  ! IF( ABS(1._SP+Z)<sqeps ) 'CLNREL : ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR -1'
  !
  rho = ABS(Z)
  IF( rho>0.375_SP ) THEN
    CLNREL = LOG(1._SP+Z)
  ELSE
    x = REAL(Z)
    CLNREL = CMPLX(0.5_SP*ALNREL(2._SP*x+rho**2),CARG(1._SP+Z),SP)
  END IF

  RETURN
END FUNCTION CLNREL