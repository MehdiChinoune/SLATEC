!** CBETA
COMPLEX FUNCTION CBETA(A,B)
  !>
  !  Compute the complete Beta function.
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
  ! CBETA computes the complete beta function of complex parameters A
  ! and B.
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
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : XERMSG
  COMPLEX :: A, B
  REAL :: xmin
  REAL, SAVE :: xmax = 0.0
  !* FIRST EXECUTABLE STATEMENT  CBETA
  IF ( xmax==0.0 ) THEN
    CALL GAMLIM(xmin,xmax)
  END IF
  !
  IF ( REAL(A)<=0.0.OR.REAL(B)<=0.0 ) CALL XERMSG('CBETA',&
    'REAL PART OF BOTH ARGUMENTS MUST BE GT 0',1,2)
  !
  IF ( REAL(A)+REAL(B)<xmax ) THEN
    CBETA = CGAMMA(A)*(CGAMMA(B)/CGAMMA(A+B))
    RETURN
  END IF
  !
  CBETA = EXP(CLBETA(A,B))
  !
END FUNCTION CBETA
