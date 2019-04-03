!** R9LN2R
REAL FUNCTION R9LN2R(X)
  IMPLICIT NONE
  !>
  !***
  !  Evaluate LOG(1+X) from second order relative accuracy so
  !            that LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X).
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

  REAL CSEVL, eps, R1MACH, sqeps, txbig, txmax, X
  INTEGER INITS
  INTEGER, SAVE :: ntln21, ntln22
  REAL, SAVE :: xmin, xbig, xmax
  REAL, PARAMETER :: ln21cs(26) = [ .18111962513478810E0,-.15627123192872463E0, &
    .028676305361557275E0, -.005558699655948139E0, .001117897665229983E0, &
    -.000230805089823279E0, .000048598853341100E0,-.000010390127388903E0, &
    .000002248456370739E0, -.000000491405927392E0, .000000108282565070E0, &
    -.000000024025872763E0, .000000005362460047E0,-.000000001202995136E0, &
    .000000000271078892E0, -.000000000061323562E0, .000000000013920858E0, &
    -.000000000003169930E0, .000000000000723837E0,-.000000000000165700E0, &
    .000000000000038018E0, -.000000000000008741E0, .000000000000002013E0, &
    -.000000000000000464E0, .000000000000000107E0, -.000000000000000024E0 ]
  REAL, PARAMETER :: ln22cs(20) = [ -.22242532535020461E0,-.061047100108078624E0, &
    .007427235009750394E0, -.000933501826163697E0, .000120049907687260E0, &
    -.000015704722952820E0, .000002081874781051E0,-.000000278919557764E0, &
    .000000037693558237E0, -.000000005130902896E0, .000000000702714117E0, &
    -.000000000096748595E0, .000000000013381046E0,-.000000000001858102E0, &
    .000000000000258929E0, -.000000000000036195E0, .000000000000005074E0, &
    -.000000000000000713E0, .000000000000000100E0,-.000000000000000014E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  R9LN2R
  IF ( first ) THEN
    eps = R1MACH(3)
    ntln21 = INITS(ln21cs,26,0.1*eps)
    ntln22 = INITS(ln22cs,20,0.1*eps)
    !
    xmin = -1.0 + SQRT(R1MACH(4))
    sqeps = SQRT(eps)
    txmax = 6.0/sqeps
    xmax = txmax - (eps*txmax**2-2.0*LOG(txmax))/(2.0*eps*txmax)
    txbig = 4.0/SQRT(sqeps)
    xbig = txbig - (sqeps*txbig**2-2.0*LOG(txbig))/(2.*sqeps*txbig)
    first = .FALSE.
  ENDIF
  !
  IF ( X<(-0.625).OR.X>0.8125 ) THEN
    !
    IF ( X<xmin ) CALL XERMSG('SLATEC','R9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1',1,1)
    IF ( X>xmax ) CALL XERMSG('SLATEC','R9LN2R',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',3,2)
    IF ( X>xbig ) CALL XERMSG('SLATEC','R9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG',2,1)
    !
    R9LN2R = (LOG(1.0+X)-X*(1.0-0.5*X))/X**3
    RETURN
  ENDIF
  !
  IF ( X<0.0 ) THEN
    R9LN2R = 0.375 + CSEVL(16.*X/5.+1.0,ln21cs,ntln21)
  ELSE
    R9LN2R = 0.375 + CSEVL(32.*X/13.-1.0,ln22cs,ntln22)
  ENDIF
  RETURN
END FUNCTION R9LN2R
