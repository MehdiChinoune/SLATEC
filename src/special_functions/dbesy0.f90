!** DBESY0
REAL(8) FUNCTION DBESY0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the Bessel function of the second kind of order
  !            zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10A1
  !***
  ! **Type:**      DOUBLE PRECISION (BESY0-S, DBESY0-D)
  !***
  ! **Keywords:**  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESY0(X) calculates the double precision Bessel function of the
  ! second kind of order zero for double precision argument X.
  !
  ! Series for BY0        on the interval  0.          to  1.60000E+01
  !                                        with weighted error   8.14E-32
  !                                         log weighted error  31.09
  !                               significant figures required  30.31
  !                                    decimal places required  31.73
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9B0MP, DBESJ0, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  INTEGER INITDS
  REAL(8) :: X, ampl, theta, y, D1MACH, DCSEVL, DBESJ0
  INTEGER, SAVE :: nty0
  REAL(8), SAVE :: xsml
  REAL(8), PARAMETER :: by0cs(19) = [ -.1127783939286557321793980546028D-1, &
    -.1283452375604203460480884531838D+0, -.1043788479979424936581762276618D+0, &
    +.2366274918396969540924159264613D-1, -.2090391647700486239196223950342D-2, &
    +.1039754539390572520999246576381D-3, -.3369747162423972096718775345037D-5, &
    +.7729384267670667158521367216371D-7, -.1324976772664259591443476068964D-8, &
    +.1764823261540452792100389363158D-10, -.1881055071580196200602823012069D-12, &
    +.1641865485366149502792237185749D-14, -.1195659438604606085745991006720D-16, &
    +.7377296297440185842494112426666D-19, -.3906843476710437330740906666666D-21, &
    +.1795503664436157949829120000000D-23, -.7229627125448010478933333333333D-26, &
    +.2571727931635168597333333333333D-28, -.8141268814163694933333333333333D-31 ]
  REAL(8), PARAMETER :: twodpi = 0.636619772367581343075535053490057D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBESY0
  IF ( first ) THEN
    nty0 = INITDS(by0cs,19,0.1*REAL(D1MACH(3)))
    xsml = SQRT(4.0D0*D1MACH(3))
    first = .FALSE.
  END IF
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBESY0','X IS ZERO OR NEGATIVE',1,2)
  IF ( X>4.0D0 ) THEN
    !
    CALL D9B0MP(X,ampl,theta)
    DBESY0 = ampl*SIN(theta)
    RETURN
  END IF
  !
  y = 0.D0
  IF ( X>xsml ) y = X*X
  DBESY0 = twodpi*LOG(0.5D0*X)*DBESJ0(X) + .375D0 + DCSEVL(.125D0*y-1.D0,by0cs,nty0)
  RETURN
END FUNCTION DBESY0
