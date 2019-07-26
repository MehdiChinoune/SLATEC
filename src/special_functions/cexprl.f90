!** CEXPRL
COMPLEX(SP) ELEMENTAL FUNCTION CEXPRL(Z)
  !> Calculate the relative error exponential (EXP(X)-1)/X.
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
  USE service, ONLY : eps_2_sp
  !
  COMPLEX(SP), INTENT(IN) :: Z
  !
  INTEGER :: i
  REAL(SP) :: r
  REAL(SP), PARAMETER :: alneps = LOG(eps_2_sp), xn = 3.72_SP - 0.3_SP*alneps, &
    xln = LOG((xn+1._SP)/1.36_SP), rbnd = eps_2_sp
  INTEGER, PARAMETER :: nterms = INT( xn - (xn*xln+alneps)/(xln+1.36_SP) + 1.5_SP )
  !* FIRST EXECUTABLE STATEMENT  CEXPRL
  !
  r = ABS(Z)
  IF( r>0.5_SP ) THEN
    CEXPRL = (EXP(Z)-1._SP)/Z
  ELSEIF( r<rbnd ) THEN
    CEXPRL = (1._SP,0._SP)
  ELSE
    CEXPRL = (0._SP,0._SP)
    DO i = 1, nterms
      CEXPRL = 1._SP + CEXPRL*Z/(nterms+2-i)
    END DO
  END IF
  !
END FUNCTION CEXPRL