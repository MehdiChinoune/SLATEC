!** R9LN2R
REAL(SP) ELEMENTAL FUNCTION R9LN2R(X)
  !> Evaluate LOG(1+X) from second order relative accuracy so that
  !            LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      SINGLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so
  ! that    LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X)
  !
  ! Series for LN21       on the interval -6.25000D-01 to  0.
  !                                        with weighted error   2.49E-17
  !                                         log weighted error  16.60
  !                               significant figures required  15.87
  !                                    decimal places required  17.31
  !
  ! Series for LN22       on the interval  0.          to  8.12500D-01
  !                                        with weighted error   1.42E-17
  !                                         log weighted error  16.85
  !                               significant figures required  15.95
  !                                    decimal places required  17.50
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : eps_2_sp, eps_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  INTEGER, PARAMETER :: ntln21 = 12, ntln22 = 9
  REAL(SP), PARAMETER :: eps = eps_2_sp, sqeps = SQRT(eps), txbig = 4._SP/SQRT(sqeps), &
    xbig = txbig - (sqeps*txbig**2-2._SP*LOG(txbig))/(2._SP*sqeps*txbig), &
    txmax = 6._SP/sqeps, xmin = -1._SP + SQRT(eps_sp), &
    xmax = txmax - (eps*txmax**2-2._SP*LOG(txmax))/(2._SP*eps*txmax)
  REAL(SP), PARAMETER :: ln21cs(26) = [ .18111962513478810_SP,-.15627123192872463_SP, &
    .028676305361557275_SP, -.005558699655948139_SP, .001117897665229983_SP, &
    -.000230805089823279_SP, .000048598853341100_SP,-.000010390127388903_SP, &
    .000002248456370739_SP, -.000000491405927392_SP, .000000108282565070_SP, &
    -.000000024025872763_SP, .000000005362460047_SP,-.000000001202995136_SP, &
    .000000000271078892_SP, -.000000000061323562_SP, .000000000013920858_SP, &
    -.000000000003169930_SP, .000000000000723837_SP,-.000000000000165700_SP, &
    .000000000000038018_SP, -.000000000000008741_SP, .000000000000002013_SP, &
    -.000000000000000464_SP, .000000000000000107_SP, -.000000000000000024_SP ]
  REAL(SP), PARAMETER :: ln22cs(20) = [ -.22242532535020461_SP,-.061047100108078624_SP, &
    .007427235009750394_SP, -.000933501826163697_SP, .000120049907687260_SP, &
    -.000015704722952820_SP, .000002081874781051_SP,-.000000278919557764_SP, &
    .000000037693558237_SP, -.000000005130902896_SP, .000000000702714117_SP, &
    -.000000000096748595_SP, .000000000013381046_SP,-.000000000001858102_SP, &
    .000000000000258929_SP, -.000000000000036195_SP, .000000000000005074_SP, &
    -.000000000000000713_SP, .000000000000000100_SP,-.000000000000000014_SP ]
  !* FIRST EXECUTABLE STATEMENT  R9LN2R
  ! ntln21 = INITS(ln21cs,0.1_SP*eps)
  ! ntln22 = INITS(ln22cs,0.1_SP*eps)
  !
  IF( X>xmax ) THEN
    ERROR STOP 'R9LN2R : NO PRECISION IN ANSWER BECAUSE X IS TOO BIG'
  ELSEIF( X<(-0.625) .OR. X>0.8125 ) THEN
    ! IF( X<xmin ) 'R9LN2R : ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1'
    ! IF( X>xbig ) 'R9LN2R : ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'
    R9LN2R = (LOG(1._SP+X)-X*(1._SP-0.5_SP*X))/X**3
  ELSEIF( X<0._SP ) THEN
    R9LN2R = 0.375_SP + CSEVL(16._SP*X/5._SP+1._SP,ln21cs(1:ntln21))
  ELSE
    R9LN2R = 0.375_SP + CSEVL(32._SP*X/13._SP-1._SP,ln22cs(1:ntln22))
  END IF

  RETURN
END FUNCTION R9LN2R