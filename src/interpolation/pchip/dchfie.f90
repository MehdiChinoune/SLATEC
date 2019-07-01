!** DCHFIE
REAL(DP) ELEMENTAL FUNCTION DCHFIE(X1,X2,F1,F2,D1,D2,A,B)
  !> Evaluates integral of a single cubic for DPCHIA
  !***
  ! **Library:**   SLATEC (PCHIP)
  !***
  ! **Type:**      DOUBLE PRECISION (CHFIE-S, DCHFIE-D)
  !***
  ! **Author:**  Fritsch, F. N., (LLNL)
  !***
  ! **Description:**
  !
  !          DCHFIE:  Cubic Hermite Function Integral Evaluator.
  !
  !     Called by  DPCHIA  to evaluate the integral of a single cubic (in
  !     Hermite form) over an arbitrary interval (A,B).
  !
  ! ----------------------------------------------------------------------
  !
  !  Calling sequence:
  !
  !        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
  !        DOUBLE PRECISION  VALUE, DCHFIE
  !
  !        VALUE = DCHFIE (X1, X2, F1, F2, D1, D2, A, B)
  !
  !   Parameters:
  !
  !     VALUE -- (output) value of the requested integral.
  !
  !     X1,X2 -- (input) endpoints if interval of definition of cubic.
  !
  !     F1,F2 -- (input) function values at the ends of the interval.
  !
  !     D1,D2 -- (input) derivative values at the ends of the interval.
  !
  !     A,B -- (input) endpoints of interval of integration.
  !
  !***
  ! **See also:**  DPCHIA
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820730  DATE WRITTEN
  !   820805  Converted to SLATEC library version.
  !   870707  Corrected subroutine name from DCHIV to DCHFIV.
  !   870813  Minor cosmetic changes.
  !   890411  1. Added SAVE statements (Vers. 3.2).
  !           2. Added SIX to DOUBLE PRECISION declaration.
  !   890411  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated AUTHOR section in prologue.  (WRB)
  !   930503  Corrected to set VALUE=0 when IERR/=0.  (FNF)
  !   930504  Eliminated IERR and changed name DCHFIV to DCHFIE.  (FNF)

  !
  !  Programming notes:
  !  1. There is no error return from this routine because zero is
  !     indeed the mathematically correct answer when X1=X2 .
  !**End
  !
  !  DECLARE ARGUMENTS.
  !
  REAL(DP), INTENT(IN) :: X1, X2, F1, F2, D1, D2, A, B
  !
  !  DECLARE LOCAL VARIABLES.
  !
  REAL(DP) :: dterm, fterm, h, phia1, phia2, phib1, phib2, psia1, psia2, psib1, &
    psib2, ta1, ta2, tb1, tb2, ua1, ua2, ub1, ub2
  !
  !  INITIALIZE.
  !
  REAL(DP), PARAMETER :: half = .5_DP, two = 2._DP, three = 3._DP, four = 4._DP, six = 6._DP
  !
  !  VALIDITY CHECK INPUT.
  !
  !* FIRST EXECUTABLE STATEMENT  DCHFIE
  IF( X1==X2 ) THEN
    DCHFIE = 0
  ELSE
    h = X2 - X1
    ta1 = (A-X1)/h
    ta2 = (X2-A)/h
    tb1 = (B-X1)/h
    tb2 = (X2-B)/h
    !
    ua1 = ta1**3
    phia1 = ua1*(two-ta1)
    psia1 = ua1*(three*ta1-four)
    ua2 = ta2**3
    phia2 = ua2*(two-ta2)
    psia2 = -ua2*(three*ta2-four)
    !
    ub1 = tb1**3
    phib1 = ub1*(two-tb1)
    psib1 = ub1*(three*tb1-four)
    ub2 = tb2**3
    phib2 = ub2*(two-tb2)
    psib2 = -ub2*(three*tb2-four)
    !
    fterm = F1*(phia2-phib2) + F2*(phib1-phia1)
    dterm = (D1*(psia2-psib2)+D2*(psib1-psia1))*(h/six)
    !
    DCHFIE = (half*h)*(fterm+dterm)
  END IF
  !
  !------------- LAST LINE OF DCHFIE FOLLOWS -----------------------------
END FUNCTION DCHFIE