!** CBETA
COMPLEX(SP) ELEMENTAL FUNCTION CBETA(A,B)
  !> Compute the complete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      COMPLEX (BETA-S, DBETA-D, CBETA-C)
  !***
  ! **Keywords:**  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CBETA computes the complete beta function of complex parameters A and B.
  ! Input Parameters:
  !       A   complex and the real part of A positive
  !       B   complex and the real part of B positive
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CGAMMA, CLBETA, GAMLIM, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section. (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : XERMSG
  COMPLEX(SP), INTENT(IN) :: A, B
  REAL(SP), PARAMETER :: xmax = 35.0307808_SP
  !* FIRST EXECUTABLE STATEMENT  CBETA
  !
  IF( REAL(A)<=0._SP .OR. REAL(B)<=0._SP ) THEN
    ERROR STOP 'CBETA : REAL PART OF BOTH ARGUMENTS MUST BE > 0'
  END IF
  !
  IF( REAL(A)+REAL(B)<xmax ) THEN
    CBETA = CGAMMA(A)*(CGAMMA(B)/CGAMMA(A+B))
  ELSE
    CBETA = EXP(CLBETA(A,B))
  END IF

  RETURN
END FUNCTION CBETA