!*==DEXPRL.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DEXPRL
DOUBLE PRECISION FUNCTION DEXPRL(X)
  IMPLICIT NONE
  !*--DEXPRL5
  !*** Start of declarations inserted by SPAG
  INTEGER i , nterms
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DEXPRL
  !***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4B
  !***TYPE      DOUBLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
  ! Taylor series is used.  If X is negative the reflection formula
  !         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
  ! may be used.  This reflection formula will be of use when the
  ! evaluation for small ABS(X) is done by Chebyshev series rather than
  ! Taylor series.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DEXPRL
  DOUBLE PRECISION X , absx , alneps , xbnd , xln , xn , D1MACH
  LOGICAL first
  SAVE nterms , xbnd , first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DEXPRL
  IF ( first ) THEN
    alneps = LOG(D1MACH(3))
    xn = 3.72D0 - 0.3D0*alneps
    xln = LOG((xn+1.0D0)/1.36D0)
    nterms = xn - (xn*xln+alneps)/(xln+1.36D0) + 1.5D0
    xbnd = D1MACH(3)
  ENDIF
  first = .FALSE.
  !
  absx = ABS(X)
  IF ( absx>0.5D0 ) DEXPRL = (EXP(X)-1.0D0)/X
  IF ( absx>0.5D0 ) RETURN
  !
  DEXPRL = 1.0D0
  IF ( absx<xbnd ) RETURN
  !
  DEXPRL = 0.0D0
  DO i = 1 , nterms
    DEXPRL = 1.0D0 + DEXPRL*X/(nterms+2-i)
  ENDDO
  !
END FUNCTION DEXPRL
