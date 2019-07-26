!** BETA
REAL(SP) ELEMENTAL FUNCTION BETA(A,B)
  !> Compute the complete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      SINGLE PRECISION (BETA-S, DBETA-D, CBETA-C)
  !***
  ! **Keywords:**  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BETA computes the complete beta function.
  !
  ! Input Parameters:
  !       A   real and positive
  !       B   real and positive
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALBETA, GAMLIM, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : tiny_sp
  !
  REAL(SP), INTENT(IN) :: A, B
  !
  REAL(SP), PARAMETER :: xmax = 35.0307808_SP
  REAL(SP), PARAMETER :: alnsml = LOG(tiny_sp)
  !* FIRST EXECUTABLE STATEMENT
  !
  IF( A<=0. .OR. B<=0. ) ERROR STOP 'BETA : BOTH ARGUMENTS MUST BE > 0'
  !
  IF( A+B<xmax ) THEN
    BETA = GAMMA(A)*GAMMA(B)/GAMMA(A+B)
  ELSE
    BETA = ALBETA(A,B)
    IF( BETA<alnsml ) THEN
      BETA = 0._SP
      ! CALL XERMSG('BETA : A AND/OR B SO BIG BETA UNDERFLOWS'
    ELSE
      BETA = EXP(BETA)
    END IF
  END IF
  !
END FUNCTION BETA