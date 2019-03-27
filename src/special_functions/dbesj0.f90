!** DBESJ0
REAL(8) FUNCTION DBESJ0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the Bessel function of the first kind of order
  !            zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10A1
  !***
  ! **Type:**      DOUBLE PRECISION (BESJ0-S, DBESJ0-D)
  !***
  ! **Keywords:**  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESJ0(X) calculates the double precision Bessel function of
  ! the first kind of order zero for double precision argument X.
  !
  ! Series for BJ0        on the interval  0.          to  1.60000E+01
  !                                        with weighted error   4.39E-32
  !                                         log weighted error  31.36
  !                               significant figures required  31.21
  !                                    decimal places required  32.00
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9B0MP, DCSEVL, INITDS

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER INITDS, ntj0
  REAL(8) :: X, ampl, theta, xsml, y, D1MACH, DCSEVL
  SAVE ntj0, xsml
  REAL(8), PARAMETER :: bj0cs(19) = [ +.10025416196893913701073127264074D+0, &
    -.66522300776440513177678757831124D+0, +.24898370349828131370460468726680D+0, &
    -.33252723170035769653884341503854D-1, +.23114179304694015462904924117729D-2, &
    -.99112774199508092339048519336549D-4, +.28916708643998808884733903747078D-5, &
    -.61210858663032635057818407481516D-7, +.98386507938567841324768748636415D-9, &
    -.12423551597301765145515897006836D-10, +.12654336302559045797915827210363D-12, &
    -.10619456495287244546914817512959D-14, +.74706210758024567437098915584000D-17, &
    -.44697032274412780547627007999999D-19, +.23024281584337436200523093333333D-21, &
    -.10319144794166698148522666666666D-23, +.40608178274873322700800000000000D-26, &
    -.14143836005240913919999999999999D-28, +.43910905496698880000000000000000D-31 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBESJ0
  IF ( first ) THEN
    ntj0 = INITDS(bj0cs,19,0.1*REAL(D1MACH(3)))
    xsml = SQRT(8.0D0*D1MACH(3))
    first = .FALSE.
  ENDIF
  !
  y = ABS(X)
  IF ( y>4.0D0 ) THEN
    !
    CALL D9B0MP(y,ampl,theta)
    DBESJ0 = ampl*COS(theta)
    RETURN
  ENDIF
  !
  DBESJ0 = 1.0D0
  IF ( y>xsml ) DBESJ0 = DCSEVL(.125D0*y*y-1.D0,bj0cs,ntj0)
  RETURN
END FUNCTION DBESJ0
