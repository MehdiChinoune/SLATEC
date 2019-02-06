!*==CBETA.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK CBETA
      COMPLEX FUNCTION CBETA(A,B)
      IMPLICIT NONE
!*--CBETA5
!*** Start of declarations inserted by SPAG
      REAL xmax , xmaxt , xmin
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  CBETA
!***PURPOSE  Compute the complete Beta function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7B
!***TYPE      COMPLEX (BETA-S, DBETA-D, CBETA-C)
!***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CBETA computes the complete beta function of complex parameters A
! and B.
! Input Parameters:
!       A   complex and the real part of A positive
!       B   complex and the real part of B positive
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CGAMMA, CLBETA, GAMLIM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890206  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  CBETA
      COMPLEX A , B , CGAMMA , CLBETA
      EXTERNAL CGAMMA
      SAVE xmax
      DATA xmax/0.0/
!***FIRST EXECUTABLE STATEMENT  CBETA
      IF ( xmax==0.0 ) THEN
        CALL GAMLIM(xmin,xmaxt)
        xmax = xmaxt
      ENDIF
!
      IF ( REAL(A)<=0.0.OR.REAL(B)<=0.0 ) CALL XERMSG('SLATEC','CBETA',
     &     'REAL PART OF BOTH ARGUMENTS MUST BE GT 0',1,2)
!
      IF ( REAL(A)+REAL(B)<xmax ) CBETA = CGAMMA(A)*(CGAMMA(B)/CGAMMA(A+B))
      IF ( REAL(A)+REAL(B)<xmax ) RETURN
!
      CBETA = EXP(CLBETA(A,B))
!
      END FUNCTION CBETA
