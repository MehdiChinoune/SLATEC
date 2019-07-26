!** DBESI0
REAL(DP) ELEMENTAL FUNCTION DBESI0(X)
  !> Compute the hyperbolic Bessel function of the first kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      DOUBLE PRECISION (BESI0-S, DBESI0-D)
  !***
  ! **Keywords:**  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESI0(X) calculates the double precision modified (hyperbolic)
  ! Bessel function of the first kind of order zero and double precision argument X.
  !
  ! Series for BI0        on the interval  0.          to  9.00000E+00
  !                                        with weighted error   9.51E-34
  !                                         log weighted error  33.02
  !                               significant figures required  33.31
  !                                    decimal places required  33.65
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DBSI0E, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : eps_2_dp, huge_dp
  !
  REAL(DP), INTENT(IN) :: X
  !
  REAL(DP) :: y
  INTEGER, PARAMETER :: nti0 = 11
  REAL(DP), PARAMETER :: xsml = SQRT(4.5_DP*eps_2_dp), xmax = LOG(huge_dp)
  REAL(DP), PARAMETER :: bi0cs(18) = [ -.7660547252839144951081894976243285E-1_DP, &
    +.1927337953993808269952408750881196E+1_DP, +.2282644586920301338937029292330415E+0_DP, &
    +.1304891466707290428079334210691888E-1_DP, +.4344270900816487451378682681026107E-3_DP, &
    +.9422657686001934663923171744118766E-5_DP, +.1434006289510691079962091878179957E-6_DP, &
    +.1613849069661749069915419719994611E-8_DP, +.1396650044535669699495092708142522E-10_DP, &
    +.9579451725505445344627523171893333E-13_DP, +.5333981859862502131015107744000000E-15_DP, &
    +.2458716088437470774696785919999999E-17_DP, +.9535680890248770026944341333333333E-20_DP, &
    +.3154382039721427336789333333333333E-22_DP, +.9004564101094637431466666666666666E-25_DP, &
    +.2240647369123670016000000000000000E-27_DP, +.4903034603242837333333333333333333E-30_DP, &
    +.9508172606122666666666666666666666E-33_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBESI0
  ! nti0 = INITDS(bi0cs,0.1_SP*eps_2_dp)
  !
  y = ABS(X)
  IF( y>xmax ) THEN
    ERROR STOP 'DBESI0 : ABS(X) SO BIG I0 OVERFLOWS'
  ELSEIF( y>3._DP ) THEN
    DBESI0 = EXP(y)*DBSI0E(X)
  ELSEIF( y>xsml ) THEN
    DBESI0 = 2.75_DP + DCSEVL(y*y/4.5_DP-1._DP,bi0cs(1:nti0))
  ELSE
    DBESI0 = 1._DP
  END IF

  RETURN
END FUNCTION DBESI0