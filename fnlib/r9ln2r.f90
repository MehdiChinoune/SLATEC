!DECK R9LN2R
FUNCTION R9LN2R(X)
  IMPLICIT NONE
  REAL CSEVL, eps, R1MACH, R9LN2R, sqeps, txbig, txmax, X, xbig, &
    xmax, xmin
  INTEGER INITS, ntln21, ntln22
  !***BEGIN PROLOGUE  R9LN2R
  !***SUBSIDIARY
  !***PURPOSE  Evaluate LOG(1+X) from second order relative accuracy so
  !            that LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4B
  !***TYPE      SINGLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  R9LN2R
  REAL ln21cs(26), ln22cs(20)
  LOGICAL first
  SAVE ln21cs, ln22cs, ntln21, ntln22, xmin, xbig, xmax, first
  DATA ln21cs(1)/.18111962513478810E0/
  DATA ln21cs(2)/ - .15627123192872463E0/
  DATA ln21cs(3)/.028676305361557275E0/
  DATA ln21cs(4)/ - .005558699655948139E0/
  DATA ln21cs(5)/.001117897665229983E0/
  DATA ln21cs(6)/ - .000230805089823279E0/
  DATA ln21cs(7)/.000048598853341100E0/
  DATA ln21cs(8)/ - .000010390127388903E0/
  DATA ln21cs(9)/.000002248456370739E0/
  DATA ln21cs(10)/ - .000000491405927392E0/
  DATA ln21cs(11)/.000000108282565070E0/
  DATA ln21cs(12)/ - .000000024025872763E0/
  DATA ln21cs(13)/.000000005362460047E0/
  DATA ln21cs(14)/ - .000000001202995136E0/
  DATA ln21cs(15)/.000000000271078892E0/
  DATA ln21cs(16)/ - .000000000061323562E0/
  DATA ln21cs(17)/.000000000013920858E0/
  DATA ln21cs(18)/ - .000000000003169930E0/
  DATA ln21cs(19)/.000000000000723837E0/
  DATA ln21cs(20)/ - .000000000000165700E0/
  DATA ln21cs(21)/.000000000000038018E0/
  DATA ln21cs(22)/ - .000000000000008741E0/
  DATA ln21cs(23)/.000000000000002013E0/
  DATA ln21cs(24)/ - .000000000000000464E0/
  DATA ln21cs(25)/.000000000000000107E0/
  DATA ln21cs(26)/ - .000000000000000024E0/
  DATA ln22cs(1)/ - .22242532535020461E0/
  DATA ln22cs(2)/ - .061047100108078624E0/
  DATA ln22cs(3)/.007427235009750394E0/
  DATA ln22cs(4)/ - .000933501826163697E0/
  DATA ln22cs(5)/.000120049907687260E0/
  DATA ln22cs(6)/ - .000015704722952820E0/
  DATA ln22cs(7)/.000002081874781051E0/
  DATA ln22cs(8)/ - .000000278919557764E0/
  DATA ln22cs(9)/.000000037693558237E0/
  DATA ln22cs(10)/ - .000000005130902896E0/
  DATA ln22cs(11)/.000000000702714117E0/
  DATA ln22cs(12)/ - .000000000096748595E0/
  DATA ln22cs(13)/.000000000013381046E0/
  DATA ln22cs(14)/ - .000000000001858102E0/
  DATA ln22cs(15)/.000000000000258929E0/
  DATA ln22cs(16)/ - .000000000000036195E0/
  DATA ln22cs(17)/.000000000000005074E0/
  DATA ln22cs(18)/ - .000000000000000713E0/
  DATA ln22cs(19)/.000000000000000100E0/
  DATA ln22cs(20)/ - .000000000000000014E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  R9LN2R
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
  ENDIF
  first = .FALSE.
  !
  IF ( X<(-0.625).OR.X>0.8125 ) THEN
    !
    IF ( X<xmin ) CALL XERMSG('SLATEC','R9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1'&
      ,1,1)
    IF ( X>xmax ) CALL XERMSG('SLATEC','R9LN2R',&
      'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG',&
      3,2)
    IF ( X>xbig ) CALL XERMSG('SLATEC','R9LN2R',&
      'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG'&
      ,2,1)
    !
    R9LN2R = (LOG(1.0+X)-X*(1.0-0.5*X))/X**3
    RETURN
  ENDIF
  !
  IF ( X<0.0 ) R9LN2R = 0.375 + CSEVL(16.*X/5.+1.0,ln21cs,ntln21)
  IF ( X>=0.0 ) R9LN2R = 0.375 + CSEVL(32.*X/13.-1.0,ln22cs,ntln22)
  RETURN
END FUNCTION R9LN2R
