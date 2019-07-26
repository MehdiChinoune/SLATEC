!** DBESK1
REAL(DP) ELEMENTAL FUNCTION DBESK1(X)
  !> Compute the modified (hyperbolic) Bessel function of the third kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      DOUBLE PRECISION (BESK1-S, DBESK1-D)
  !***
  ! **Keywords:**  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESK1(X) calculates the double precision modified (hyperbolic)
  ! Bessel function of the third kind of order one for double precision
  ! argument X.  The argument must be large enough that the result does
  ! not overflow and small enough that the result does not underflow.
  !
  ! Series for BK1        on the interval  0.          to  4.00000E+00
  !                                        with weighted error   9.16E-32
  !                                         log weighted error  31.04
  !                               significant figures required  30.61
  !                                    decimal places required  31.64
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : tiny_dp, huge_dp, eps_2_dp
  !
  REAL(DP), INTENT(IN) :: X
  !
  REAL(DP) :: y
  INTEGER, PARAMETER :: ntk1 = 10
  REAL(DP), PARAMETER :: xmin = EXP(MAX(LOG(tiny_dp),-LOG(huge_dp))+0.01_DP), &
    xsml = SQRT(4._DP*eps_2_dp), xmaxt = -LOG(tiny_dp), &
    xmax = xmaxt - 0.5_DP*xmaxt*LOG(xmaxt)/(xmaxt+0.5_DP)
  REAL(DP), PARAMETER :: bk1cs(16) = [ +.25300227338947770532531120868533E-1_DP, &
    -.35315596077654487566723831691801E+0_DP,-.12261118082265714823479067930042E+0_DP, &
    -.69757238596398643501812920296083E-2_DP,-.17302889575130520630176507368979E-3_DP, &
    -.24334061415659682349600735030164E-5_DP,-.22133876307347258558315252545126E-7_DP, &
    -.14114883926335277610958330212608E-9_DP,-.66669016941993290060853751264373E-12_DP, &
    -.24274498505193659339263196864853E-14_DP,-.70238634793862875971783797120000E-17_DP, &
    -.16543275155100994675491029333333E-19_DP,-.32338347459944491991893333333333E-22_DP, &
    -.53312750529265274999466666666666E-25_DP,-.75130407162157226666666666666666E-28_DP, &
    -.91550857176541866666666666666666E-31_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBESK1
  ! ntk1 = INITDS(bk1cs,0.1_SP*eps_2_dp)
  !
  IF( X<=0._DP ) THEN
    ERROR STOP 'DBESK1 : X IS ZERO OR NEGATIVE'
  ELSEIF( X<xmin ) THEN
    ERROR STOP 'DBESK1 : X SO SMALL K1 OVERFLOWS'
  ELSEIF( X<=2._DP ) THEN
    y = 0._DP
    IF( X>xsml ) y = X*X
    DBESK1 = LOG(0.5_DP*X)*DBESI1(X) + (0.75_DP+DCSEVL(.5_DP*y-1._DP,bk1cs(1:ntk1)))/X
  ELSEIF( X<=xmax ) THEN
    DBESK1 = EXP(-X)*DBSK1E(X)
  ELSE
    DBESK1 = 0._DP
    ! IF( X>xmax ) CALL XERMSG('DBESK1','X SO BIG K1 UNDERFLOWS',1,1)
  END IF

  RETURN
END FUNCTION DBESK1