!** BI
REAL FUNCTION BI(X)
  !>
  !***
  !  Evaluate the Bairy function (the Airy function of the
  !            second kind).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      SINGLE PRECISION (BI-S, DBI-D)
  !***
  ! **Keywords:**  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BI(X) calculates the Airy function of the second kind for real
  ! argument X.
  !
  ! Series for BIF        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   1.88E-19
  !                                         log weighted error  18.72
  !                               significant figures required  17.74
  !                                    decimal places required  19.20
  !
  ! Series for BIG        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   2.61E-17
  !                                         log weighted error  16.58
  !                               significant figures required  15.17
  !                                    decimal places required  17.03
  !
  ! Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00
  !                                        with weighted error   1.11E-17
  !                                         log weighted error  16.95
  !                        approx significant figures required  16.5
  !                                    decimal places required  17.45
  !
  ! Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00
  !                                        with weighted error   1.19E-18
  !                                         log weighted error  17.92
  !                        approx significant figures required  17.2
  !                                    decimal places required  18.42
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL eta, theta, X, xm, z
  INTEGER, SAVE :: nbif, nbig, nbif2, nbig2
  REAL, SAVE :: x3sml, xmax
  REAL, PARAMETER :: bifcs(9) = [ -.01673021647198664948E0, .1025233583424944561E0, &
    .00170830925073815165E0, .00001186254546774468E0, .00000004493290701779E0, &
    .00000000010698207143E0, .00000000000017480643E0, .00000000000000020810E0, &
    .00000000000000000018E0 ]
  REAL, PARAMETER :: bigcs(8) = [ .02246622324857452E0, .03736477545301955E0, &
    .00044476218957212E0, .00000247080756363E0, .00000000791913533E0, &
    .00000000001649807E0, .00000000000002411E0, .00000000000000002E0 ]
  REAL, PARAMETER :: bif2cs(10) = [ 0.09984572693816041E0, .478624977863005538E0, &
    .0251552119604330118E0, .0005820693885232645E0, .0000074997659644377E0, &
    .0000000613460287034E0, .0000000003462753885E0, .0000000000014288910E0, &
    .0000000000000044962E0, .0000000000000000111E0 ]
  REAL, PARAMETER :: big2cs(10) = [ .033305662145514340E0, .161309215123197068E0, &
    .0063190073096134286E0, .0001187904568162517E0, .0000013045345886200E0, &
    .0000000093741259955E0, .0000000000474580188E0, .0000000000001783107E0, &
    .0000000000000005167E0, .0000000000000000011E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BI
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    nbif = INITS(bifcs,9,eta)
    nbig = INITS(bigcs,8,eta)
    nbif2 = INITS(bif2cs,10,eta)
    nbig2 = INITS(big2cs,10,eta)
    !
    x3sml = eta**0.3333
    xmax = (1.5*LOG(R1MACH(2)))**0.6666
    first = .FALSE.
  END IF
  !
  IF ( X<(-1.0) ) THEN
    CALL R9AIMP(X,xm,theta)
    BI = xm*SIN(theta)
    RETURN
    !
  ELSEIF ( X<=1.0 ) THEN
    z = 0.0
    IF ( ABS(X)>x3sml ) z = X**3
    BI = 0.625 + CSEVL(z,bifcs,nbif) + X*(0.4375+CSEVL(z,bigcs,nbig))
    RETURN
    !
  ELSEIF ( X>2.0 ) THEN
    !
    IF ( X>xmax ) CALL XERMSG('SLATEC','BI','X SO BIG THAT BI OVERFLOWS',1,2)
    !
    BI = BIE(X)*EXP(2.0*X*SQRT(X)/3.0)
    RETURN
  END IF
  z = (2.0*X**3-9.0)/7.0
  BI = 1.125 + CSEVL(z,bif2cs,nbif2) + X*(0.625+CSEVL(z,big2cs,nbig2))
  RETURN
END FUNCTION BI
