!** BI
REAL(SP) ELEMENTAL FUNCTION BI(X)
  !> Evaluate the Bairy function (the Airy function of the second kind).
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
  ! BI(X) calculates the Airy function of the second kind for real argument X.
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
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: theta, xm, z
  INTEGER, PARAMETER :: nbif = 5, nbig = 5, nbif2 = 6, nbig2 = 6
  REAL(SP), PARAMETER :: eta = 0.1_SP*R1MACH(3), x3sml = eta**0.3333_SP, &
    xmax = (1.5_SP*LOG(R1MACH(2)))**0.6666_SP
  REAL(SP), PARAMETER :: bifcs(9) = [ -.01673021647198664948_SP, .1025233583424944561_SP, &
    .00170830925073815165_SP, .00001186254546774468_SP, .00000004493290701779_SP, &
    .00000000010698207143_SP, .00000000000017480643_SP, .00000000000000020810_SP, &
    .00000000000000000018_SP ]
  REAL(SP), PARAMETER :: bigcs(8) = [ .02246622324857452_SP, .03736477545301955_SP, &
    .00044476218957212_SP, .00000247080756363_SP, .00000000791913533_SP, &
    .00000000001649807_SP, .00000000000002411_SP, .00000000000000002_SP ]
  REAL(SP), PARAMETER :: bif2cs(10) = [ 0.09984572693816041_SP, .478624977863005538_SP, &
    .0251552119604330118_SP, .0005820693885232645_SP, .0000074997659644377_SP, &
    .0000000613460287034_SP, .0000000003462753885_SP, .0000000000014288910_SP, &
    .0000000000000044962_SP, .0000000000000000111_SP ]
  REAL(SP), PARAMETER :: big2cs(10) = [ .033305662145514340_SP, .161309215123197068_SP, &
    .0063190073096134286_SP, .0001187904568162517_SP, .0000013045345886200_SP, &
    .0000000093741259955_SP, .0000000000474580188_SP, .0000000000001783107_SP, &
    .0000000000000005167_SP, .0000000000000000011_SP ]
  !* FIRST EXECUTABLE STATEMENT  BI
  ! nbif = INITS(bifcs,eta)
  ! nbig = INITS(bigcs,eta)
  ! nbif2 = INITS(bif2cs,eta)
  ! nbig2 = INITS(big2cs,eta)
  !
  IF( X<(-1._SP) ) THEN
    CALL R9AIMP(X,xm,theta)
    BI = xm*SIN(theta)
  ELSEIF( X<=1._SP ) THEN
    z = 0._SP
    IF( ABS(X)>x3sml ) z = X**3
    BI = 0.625_SP + CSEVL(z,bifcs(1:nbif)) + X*(0.4375_SP+CSEVL(z,bigcs(1:nbig)))
  ELSEIF( x<=2._SP ) THEN
    z = (2._SP*X**3-9._SP)/7._SP
    BI = 1.125_SP + CSEVL(z,bif2cs(1:nbif2)) + X*(0.625_SP+CSEVL(z,big2cs(1:nbig2)))
  ELSEIF( X<=xmax ) THEN
    BI = BIE(X)*EXP(2._SP*X*SQRT(X)/3._SP)
  ELSE
    ERROR STOP 'BI : X SO BIG THAT BI OVERFLOWS'
  END IF

  RETURN
END FUNCTION BI