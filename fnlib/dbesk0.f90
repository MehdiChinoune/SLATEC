!*==DBESK0.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBESK0
REAL(8) FUNCTION DBESK0(X)
  IMPLICIT NONE
  !*--DBESK05
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , ntk0
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBESK0
  !***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
  !            third kind of order zero.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B1
  !***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D)
  !***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBESK0(X) calculates the double precision modified (hyperbolic)
  ! Bessel function of the third kind of order zero for double
  ! precision argument X.  The argument must be greater than zero
  ! but not so large that the result underflows.
  !
  ! Series for BK0        on the interval  0.          to  4.00000E+00
  !                                        with weighted error   3.08E-33
  !                                         log weighted error  32.51
  !                               significant figures required  32.05
  !                                    decimal places required  33.11
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBESK0
  REAL(8) :: X , bk0cs(16) , xmax , xmaxt , xsml , y , D1MACH , &
    DCSEVL , DBESI0 , DBSK0E
  LOGICAL first
  SAVE bk0cs , ntk0 , xsml , xmax , first
  DATA bk0cs(1)/ - .353273932339027687201140060063153D-1/
  DATA bk0cs(2)/ + .344289899924628486886344927529213D+0/
  DATA bk0cs(3)/ + .359799365153615016265721303687231D-1/
  DATA bk0cs(4)/ + .126461541144692592338479508673447D-2/
  DATA bk0cs(5)/ + .228621210311945178608269830297585D-4/
  DATA bk0cs(6)/ + .253479107902614945730790013428354D-6/
  DATA bk0cs(7)/ + .190451637722020885897214059381366D-8/
  DATA bk0cs(8)/ + .103496952576336245851008317853089D-10/
  DATA bk0cs(9)/ + .425981614279108257652445327170133D-13/
  DATA bk0cs(10)/ + .137446543588075089694238325440000D-15/
  DATA bk0cs(11)/ + .357089652850837359099688597333333D-18/
  DATA bk0cs(12)/ + .763164366011643737667498666666666D-21/
  DATA bk0cs(13)/ + .136542498844078185908053333333333D-23/
  DATA bk0cs(14)/ + .207527526690666808319999999999999D-26/
  DATA bk0cs(15)/ + .271281421807298560000000000000000D-29/
  DATA bk0cs(16)/ + .308259388791466666666666666666666D-32/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBESK0
  IF ( first ) THEN
    ntk0 = INITDS(bk0cs,16,0.1*REAL(D1MACH(3)))
    xsml = SQRT(4.0D0*D1MACH(3))
    xmaxt = -LOG(D1MACH(1))
    xmax = xmaxt - 0.5D0*xmaxt*LOG(xmaxt)/(xmaxt+0.5D0)
  ENDIF
  first = .FALSE.
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBESK0','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2.0D0 ) THEN
    !
    DBESK0 = 0.D0
    IF ( X>xmax ) CALL XERMSG('SLATEC','DBESK0','X SO BIG K0 UNDERFLOWS',1,&
      1)
    IF ( X>xmax ) RETURN
    !
    DBESK0 = EXP(-X)*DBSK0E(X)
    GOTO 99999
  ENDIF
  !
  y = 0.D0
  IF ( X>xsml ) y = X*X
  DBESK0 = -LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + DCSEVL(.5D0*y-1.D0,bk0cs,ntk0)
  RETURN
  !
  99999 END FUNCTION DBESK0
