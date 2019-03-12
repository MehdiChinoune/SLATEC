!DECK DBETA
REAL(8) FUNCTION DBETA(A,B)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DBETA
  !***PURPOSE  Compute the complete Beta function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7B
  !***TYPE      DOUBLE PRECISION (BETA-S, DBETA-D, CBETA-C)
  !***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBETA(A,B) calculates the double precision complete beta function
  ! for double precision arguments A and B.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DGAMLM, DGAMMA, DLBETA, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  DBETA
  REAL(8) :: A, B, alnsml, xmax, xmin, DLBETA, DGAMMA, D1MACH
  LOGICAL first
  EXTERNAL DGAMMA
  SAVE xmax, alnsml, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBETA
  IF ( first ) THEN
    CALL DGAMLM(xmin,xmax)
    alnsml = LOG(D1MACH(1))
  ENDIF
  first = .FALSE.
  !
  IF ( A<=0.D0.OR.B<=0.D0 )&
    CALL XERMSG('SLATEC','DBETA','BOTH ARGUMENTS MUST BE GT 0',2,2)
  !
  IF ( A+B<xmax ) THEN
    DBETA = DGAMMA(A)*DGAMMA(B)/DGAMMA(A+B)
    RETURN
  ENDIF
  !
  DBETA = DLBETA(A,B)
  IF ( DBETA<alnsml ) THEN
    !
    DBETA = 0.D0
    CALL XERMSG('SLATEC','DBETA','A AND/OR B SO BIG BETA UNDERFLOWS',1,1)
    RETURN
  ENDIF
  DBETA = EXP(DBETA)
  RETURN
END FUNCTION DBETA
