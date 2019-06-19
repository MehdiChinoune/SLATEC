!** ALNREL
REAL(SP) FUNCTION ALNREL(X)
  !> Evaluate ln(1+X) accurate in the sense of relative error.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      SINGLE PRECISION (ALNREL-S, DLNREL-D, CLNREL-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ALNREL(X) evaluates ln(1+X) accurately in the sense of relative
  ! error when X is very small.  This routine must be used to
  ! maintain relative error accuracy whenever X is small and
  ! accurately known.
  !
  ! Series for ALNR       on the interval -3.75000D-01 to  3.75000D-01
  !                                        with weighted error   1.93E-17
  !                                         log weighted error  16.72
  !                               significant figures required  16.44
  !                                    decimal places required  17.40
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL(SP) :: X
  INTEGER, SAVE :: nlnrel
  REAL(SP), PARAMETER :: xmin = -1.0 + SQRT(R1MACH(4))
  REAL(SP), PARAMETER :: alnrcs(23) = [ 1.0378693562743770E0, -.13364301504908918E0, &
    .019408249135520563E0, -.003010755112753577E0, .000486946147971548E0, &
    -.000081054881893175E0, .000013778847799559E0,-.000002380221089435E0, &
    .000000416404162138E0, -.000000073595828378E0, .000000013117611876E0, &
    -.000000002354670931E0, .000000000425227732E0,-.000000000077190894E0, &
    .000000000014075746E0, -.000000000002576907E0, .000000000000473424E0, &
    -.000000000000087249E0, .000000000000016124E0,-.000000000000002987E0, &
    .000000000000000554E0, -.000000000000000103E0, .000000000000000019E0 ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  ALNREL
  IF( first ) THEN
    nlnrel = INITS(alnrcs,23,0.1*R1MACH(3))
  END IF
  !
  IF( X<=(-1.0) ) CALL XERMSG('ALNREL','X IS LE -1',2,2)
  IF( X<xmin ) CALL XERMSG('ALNREL',&
    'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1',1,1)
  !
  IF( ABS(X)<=0.375 ) THEN
    ALNREL = X*(1.-X*CSEVL(X/.375,alnrcs,nlnrel))
  ELSE
    ALNREL = LOG(1.0+X)
  END IF
  !
END FUNCTION ALNREL
