!** DEXPRL
REAL(DP) FUNCTION DEXPRL(X)
  !> Calculate the relative error exponential (EXP(X)-1)/X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      DOUBLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
  ! Taylor series is used.  If X is negative the reflection formula
  !         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
  ! may be used.  This reflection formula will be of use when the
  ! evaluation for small ABS(X) is done by Chebyshev series rather than
  ! Taylor series.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  REAL(DP) :: X
  INTEGER :: i
  REAL(DP) :: absx
  REAL(DP), PARAMETER :: alneps = LOG(D1MACH(3)), xn = 3.72_DP - 0.3_DP*alneps, &
    xln = LOG((xn+1._DP)/1.36_DP), xbnd = D1MACH(3)
  INTEGER, PARAMETER :: nterms = INT( xn - (xn*xln+alneps)/(xln+1.36_DP) + 1.5_DP )
  !* FIRST EXECUTABLE STATEMENT  DEXPRL
  !
  absx = ABS(X)
  IF( absx>0.5_DP ) DEXPRL = (EXP(X)-1._DP)/X
  IF( absx>0.5_DP ) RETURN
  !
  DEXPRL = 1._DP
  IF( absx<xbnd ) RETURN
  !
  DEXPRL = 0._DP
  DO i = 1, nterms
    DEXPRL = 1._DP + DEXPRL*X/(nterms+2-i)
  END DO
  !
END FUNCTION DEXPRL
