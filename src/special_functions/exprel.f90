!** EXPREL
REAL FUNCTION EXPREL(X)
  IMPLICIT NONE
  !>
  !***
  !  Calculate the relative error exponential (EXP(X)-1)/X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      SINGLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
  ! Taylor series is used.  If X is negative, the reflection formula
  !         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
  ! may be used.  This reflection formula will be of use when the
  ! evaluation for small ABS(X) is done by Chebyshev series rather than
  ! Taylor series.  EXPREL and X are single precision.
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
  
  REAL absx, alneps, R1MACH, X, xbnd, xln, xn
  INTEGER i, nterms
  LOGICAL first
  SAVE nterms, xbnd, first
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  EXPREL
  IF ( first ) THEN
    alneps = LOG(R1MACH(3))
    xn = 3.72 - 0.3*alneps
    xln = LOG((xn+1.0)/1.36)
    nterms = INT( xn - (xn*xln+alneps)/(xln+1.36) + 1.5 )
    xbnd = R1MACH(3)
  ENDIF
  first = .FALSE.
  !
  absx = ABS(X)
  IF ( absx>0.5 ) EXPREL = (EXP(X)-1.0)/X
  IF ( absx>0.5 ) RETURN
  !
  EXPREL = 1.0
  IF ( absx<xbnd ) RETURN
  !
  EXPREL = 0.0
  DO i = 1, nterms
    EXPREL = 1.0 + EXPREL*X/(nterms+2-i)
  ENDDO
  !
END FUNCTION EXPREL
