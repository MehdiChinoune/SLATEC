!** DBESK1
REAL(8) FUNCTION DBESK1(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the modified (hyperbolic) Bessel function of the
  !            third kind of order one.
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

  INTEGER INITDS
  REAL(8) :: X, xmaxt, y, D1MACH, DCSEVL, DBESI1, DBSK1E
  INTEGER, SAVE :: ntk1
  REAL(8), SAVE :: xmin, xsml, xmax
  REAL(8), PARAMETER :: bk1cs(16) = [ +.25300227338947770532531120868533D-1, &
    -.35315596077654487566723831691801D+0,-.12261118082265714823479067930042D+0, &
    -.69757238596398643501812920296083D-2,-.17302889575130520630176507368979D-3, &
    -.24334061415659682349600735030164D-5,-.22133876307347258558315252545126D-7, &
    -.14114883926335277610958330212608D-9,-.66669016941993290060853751264373D-12, &
    -.24274498505193659339263196864853D-14,-.70238634793862875971783797120000D-17, &
    -.16543275155100994675491029333333D-19,-.32338347459944491991893333333333D-22, &
    -.53312750529265274999466666666666D-25,-.75130407162157226666666666666666D-28, &
    -.91550857176541866666666666666666D-31 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBESK1
  IF ( first ) THEN
    ntk1 = INITDS(bk1cs,16,0.1*REAL(D1MACH(3)))
    xmin = EXP(MAX(LOG(D1MACH(1)),-LOG(D1MACH(2)))+0.01D0)
    xsml = SQRT(4.0D0*D1MACH(3))
    xmaxt = -LOG(D1MACH(1))
    xmax = xmaxt - 0.5D0*xmaxt*LOG(xmaxt)/(xmaxt+0.5D0)
    first = .FALSE.
  ENDIF
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBESK1','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2.0D0 ) THEN
    !
    DBESK1 = 0.D0
    IF ( X>xmax ) CALL XERMSG('SLATEC','DBESK1','X SO BIG K1 UNDERFLOWS',1,1)
    IF ( X>xmax ) RETURN
    !
    DBESK1 = EXP(-X)*DBSK1E(X)
    RETURN
  ENDIF
  !
  IF ( X<xmin ) CALL XERMSG('SLATEC','DBESK1','X SO SMALL K1 OVERFLOWS',3,2)
  y = 0.D0
  IF ( X>xsml ) y = X*X
  DBESK1 = LOG(0.5D0*X)*DBESI1(X) + (0.75D0+DCSEVL(.5D0*y-1.D0,bk1cs,ntk1))/X
  RETURN
END FUNCTION DBESK1
