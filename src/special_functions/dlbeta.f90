!** DLBETA
REAL(8) FUNCTION DLBETA(A,B)
  !>
  !  Compute the natural logarithm of the complete Beta
  !            function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      DOUBLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C)
  !***
  ! **Keywords:**  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DLBETA(A,B) calculates the double precision natural logarithm of
  ! the complete beta function for double precision arguments
  ! A and B.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D9LGMC, DLNREL, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : XERMSG
  REAL(8) :: A, B, p, q, corr
  REAL(8), PARAMETER :: sq2pil = 0.91893853320467274178032973640562D0
  !* FIRST EXECUTABLE STATEMENT  DLBETA
  p = MIN(A,B)
  q = MAX(A,B)
  !
  IF ( p<=0.D0 ) CALL XERMSG('DLBETA',&
    'BOTH ARGUMENTS MUST BE GT ZERO',1,2)
  !
  IF ( p>=10.D0 ) THEN
    !
    ! P AND Q ARE BIG.
    !
    corr = D9LGMC(p) + D9LGMC(q) - D9LGMC(p+q)
    DLBETA = -0.5D0*LOG(q) + sq2pil + corr + (p-0.5D0)*LOG(p/(p+q))&
      + q*DLNREL(-p/(p+q))
    RETURN
  ELSEIF ( q<10.D0 ) THEN
    !
    ! P AND Q ARE SMALL.
    !
    DLBETA = LOG(GAMMA(p)*(GAMMA(q)/GAMMA(p+q)))
    RETURN
  END IF
  !
  ! P IS SMALL, BUT Q IS BIG.
  !
  corr = D9LGMC(q) - D9LGMC(p+q)
  DLBETA = LOG_GAMMA(p) + corr + p - p*LOG(p+q) + (q-0.5D0)*DLNREL(-p/(p+q))
  RETURN
END FUNCTION DLBETA
