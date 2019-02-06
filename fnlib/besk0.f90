!*==BESK0.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK BESK0
      FUNCTION BESK0(X)
      IMPLICIT NONE
!*--BESK05
!*** Start of declarations inserted by SPAG
      REAL BESI0 , BESK0 , BESK0E , bk0cs , CSEVL , R1MACH , X , xmax , xmaxt , 
     &     xsml , y
      INTEGER INITS , ntk0
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  BESK0
!***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
!            third kind of order zero.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      SINGLE PRECISION (BESK0-S, DBESK0-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BESK0(X) calculates the modified (hyperbolic) Bessel function
! of the third kind of order zero for real argument X .GT. 0.0.
!
! Series for BK0        on the interval  0.          to  4.00000D+00
!                                        with weighted error   3.57E-19
!                                         log weighted error  18.45
!                               significant figures required  17.99
!                                    decimal places required  18.97
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BESK0
      DIMENSION bk0cs(11)
      LOGICAL first
      SAVE bk0cs , ntk0 , xsml , xmax , first
      DATA bk0cs(1)/ - .03532739323390276872E0/
      DATA bk0cs(2)/.3442898999246284869E0/
      DATA bk0cs(3)/.03597993651536150163E0/
      DATA bk0cs(4)/.00126461541144692592E0/
      DATA bk0cs(5)/.00002286212103119451E0/
      DATA bk0cs(6)/.00000025347910790261E0/
      DATA bk0cs(7)/.00000000190451637722E0/
      DATA bk0cs(8)/.00000000001034969525E0/
      DATA bk0cs(9)/.00000000000004259816E0/
      DATA bk0cs(10)/.00000000000000013744E0/
      DATA bk0cs(11)/.00000000000000000035E0/
      DATA first/.TRUE./
!***FIRST EXECUTABLE STATEMENT  BESK0
      IF ( first ) THEN
        ntk0 = INITS(bk0cs,11,0.1*R1MACH(3))
        xsml = SQRT(4.0*R1MACH(3))
        xmaxt = -LOG(R1MACH(1))
        xmax = xmaxt - 0.5*xmaxt*LOG(xmaxt)/(xmaxt+0.5) - 0.01
      ENDIF
      first = .FALSE.
!
      IF ( X<=0. ) CALL XERMSG('SLATEC','BESK0','X IS ZERO OR NEGATIVE',2,2)
      IF ( X>2. ) THEN
!
        BESK0 = 0.
        IF ( X>xmax ) CALL XERMSG('SLATEC','BESK0','X SO BIG K0 UNDERFLOWS',1,1)
        IF ( X>xmax ) RETURN
!
        BESK0 = EXP(-X)*BESK0E(X)
        GOTO 99999
      ENDIF
!
      y = 0.
      IF ( X>xsml ) y = X*X
      BESK0 = -LOG(0.5*X)*BESI0(X) - .25 + CSEVL(.5*y-1.,bk0cs,ntk0)
      RETURN
!
99999 END FUNCTION BESK0
