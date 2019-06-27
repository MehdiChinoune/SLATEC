!** CLBETA
COMPLEX(SP) ELEMENTAL FUNCTION CLBETA(A,B)
  !> Compute the natural logarithm of the complete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      COMPLEX (ALBETA-S, DLBETA-D, CLBETA-C)
  !***
  ! **Keywords:**  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CLBETA computes the natural log of the complex valued complete beta
  ! function of complex parameters A and B.  This is a preliminary version
  ! which is not accurate.
  !
  ! Input Parameters:
  !       A   complex and the real part of A positive
  !       B   complex and the real part of B positive
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CLNGAM, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  COMPLEX(SP), INTENT(IN) :: A, B
  !* FIRST EXECUTABLE STATEMENT  CLBETA
  IF( REAL(A)<=0._SP .OR. REAL(B)<=0._SP ) THEN
    ERROR STOP 'CLBETA : REAL PART OF BOTH ARGUMENTS MUST BE GT 0'
  END IF
  !
  CLBETA = CLNGAM(A) + CLNGAM(B) - CLNGAM(A+B)
  !
END FUNCTION CLBETA