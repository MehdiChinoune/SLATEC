!** ALNREL
REAL(SP) ELEMENTAL FUNCTION ALNREL(X)
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
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !
  REAL(SP), INTENT(IN) :: X
  !
  INTEGER, PARAMETER :: nlnrel = 11
  REAL(SP), PARAMETER :: alnrcs(23) = [ 1.0378693562743770_SP, -.13364301504908918_SP, &
    .019408249135520563_SP, -.003010755112753577_SP, .000486946147971548_SP, &
    -.000081054881893175_SP, .000013778847799559_SP,-.000002380221089435_SP, &
    .000000416404162138_SP, -.000000073595828378_SP, .000000013117611876_SP, &
    -.000000002354670931_SP, .000000000425227732_SP,-.000000000077190894_SP, &
    .000000000014075746_SP, -.000000000002576907_SP, .000000000000473424_SP, &
    -.000000000000087249_SP, .000000000000016124_SP,-.000000000000002987_SP, &
    .000000000000000554_SP, -.000000000000000103_SP, .000000000000000019_SP ]
  !* FIRST EXECUTABLE STATEMENT  ALNREL
  ! nlnrel = INITS(alnrcs,0.1_SP*eps_2_sp)
  !
  IF( X<=(-1._SP) ) THEN
    ERROR STOP 'ALNREL : X IS <= -1'
  ! IF( X<xmin ) CALL XERMSG('ALNREL : ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1',1,1)
  ELSEIF( ABS(X)<=0.375_SP ) THEN
    ALNREL = X*(1._SP-X*CSEVL(X/.375_SP,alnrcs(1:nlnrel)))
  ELSE
    ALNREL = LOG(1._SP+X)
  END IF

  RETURN
END FUNCTION ALNREL