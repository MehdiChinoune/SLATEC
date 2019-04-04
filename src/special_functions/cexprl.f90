!** CEXPRL
COMPLEX FUNCTION CEXPRL(Z)
  IMPLICIT NONE
  !>
  !***
  !  Calculate the relative error exponential (EXP(X)-1)/X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      COMPLEX (EXPREL-S, DEXPRL-D, CEXPRL-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  (EXP(Z)-1)/Z .  For small ABS(Z), we use the Taylor
  ! series.  We could instead use the expression
  !        CEXPRL(Z) = (EXP(X)*EXP(I*Y)-1)/Z
  !                  = (X*EXPREL(X) * (1 - 2*SIN(Y/2)**2) - 2*SIN(Y/2)**2
  !                                    + I*SIN(Y)*(1+X*EXPREL(X))) / Z
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  REAL alneps, r, R1MACH, xln, xn
  INTEGER i
  COMPLEX Z
  INTEGER, SAVE :: nterms
  REAL, SAVE :: rbnd
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  CEXPRL
  IF ( first ) THEN
    alneps = LOG(R1MACH(3))
    xn = 3.72 - 0.3*alneps
    xln = LOG((xn+1.0)/1.36)
    nterms = INT( xn - (xn*xln+alneps)/(xln+1.36) + 1.5 )
    rbnd = R1MACH(3)
    first = .FALSE.
  END IF
  !
  r = ABS(Z)
  IF ( r>0.5 ) CEXPRL = (EXP(Z)-1.0)/Z
  IF ( r>0.5 ) RETURN
  !
  CEXPRL = (1.0,0.0)
  IF ( r<rbnd ) RETURN
  !
  CEXPRL = (0.0,0.0)
  DO i = 1, nterms
    CEXPRL = 1.0 + CEXPRL*Z/(nterms+2-i)
  END DO
  !
END FUNCTION CEXPRL
