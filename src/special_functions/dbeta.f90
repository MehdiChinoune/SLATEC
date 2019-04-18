!** DBETA
REAL(8) FUNCTION DBETA(A,B)
  !>
  !***
  !  Compute the complete Beta function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7B
  !***
  ! **Type:**      DOUBLE PRECISION (BETA-S, DBETA-D, CBETA-C)
  !***
  ! **Keywords:**  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBETA(A,B) calculates the double precision complete beta function
  ! for double precision arguments A and B.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DGAMLM, DLBETA, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  REAL(8) :: A, B, xmin
  REAL(8), SAVE :: xmax, alnsml
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBETA
  IF ( first ) THEN
    CALL DGAMLM(xmin,xmax)
    alnsml = LOG(D1MACH(1))
    first = .FALSE.
  END IF
  !
  IF ( A<=0.D0.OR.B<=0.D0 )&
    CALL XERMSG('SLATEC','DBETA','BOTH ARGUMENTS MUST BE GT 0',2,2)
  !
  IF ( A+B<xmax ) THEN
    DBETA = GAMMA(A)*GAMMA(B)/GAMMA(A+B)
    RETURN
  END IF
  !
  DBETA = DLBETA(A,B)
  IF ( DBETA<alnsml ) THEN
    !
    DBETA = 0.D0
    CALL XERMSG('SLATEC','DBETA','A AND/OR B SO BIG BETA UNDERFLOWS',1,1)
    RETURN
  END IF
  DBETA = EXP(DBETA)
  RETURN
END FUNCTION DBETA
