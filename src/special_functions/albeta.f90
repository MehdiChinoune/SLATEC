!** ALBETA
REAL(SP) FUNCTION ALBETA(A,B)
  !>
  !  Compute the natural logarithm of the complete Beta
  !            function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      SINGLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C)
  !***
  ! **Keywords:**  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ALBETA computes the natural log of the complete beta function.
  !
  ! Input Parameters:
  !       A   real and positive
  !       B   real and positive
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALNREL, R9LGMC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : XERMSG
  REAL(SP) A, B, corr, p, q
  REAL(SP), PARAMETER :: sq2pil = 0.91893853320467274E0
  !* FIRST EXECUTABLE STATEMENT  ALBETA
  p = MIN(A,B)
  q = MAX(A,B)
  !
  IF ( p<=0.0 ) CALL XERMSG('ALBETA',&
    'BOTH ARGUMENTS MUST BE GT ZERO',1,2)
  IF ( p>=10.0 ) THEN
    !
    ! P AND Q ARE BIG.
    !
    corr = R9LGMC(p) + R9LGMC(q) - R9LGMC(p+q)
    ALBETA = -0.5*LOG(q) + sq2pil + corr + (p-0.5)*LOG(p/(p+q))+ q*ALNREL(-p/(p+q))
    RETURN
  ELSEIF ( q<10.0 ) THEN
    !
    ! P AND Q ARE SMALL.
    !
    ALBETA = LOG(GAMMA(p)*(GAMMA(q)/GAMMA(p+q)))
    RETURN
  END IF
  !
  ! P IS SMALL, BUT Q IS BIG.
  !
  corr = R9LGMC(q) - R9LGMC(p+q)
  ALBETA = LOG_GAMMA(p) + corr + p - p*LOG(p+q) + (q-0.5)*ALNREL(-p/(p+q))
  RETURN
END FUNCTION ALBETA
