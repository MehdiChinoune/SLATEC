!** DERF
REAL(8) FUNCTION DERF(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the error function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C8A, L5A1E
  !***
  ! **Type:**      DOUBLE PRECISION (ERF-S, DERF-D)
  !***
  ! **Keywords:**  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DERF(X) calculates the double precision error function for double
  ! precision argument X.
  !
  ! Series for ERF        on the interval  0.          to  1.00000E+00
  !                                        with weighted error   1.28E-32
  !                                         log weighted error  31.89
  !                               significant figures required  31.05
  !                                    decimal places required  32.55
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, DERFC, INITDS

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable name.  (RWC, WRB)

  INTEGER nterf
  REAL(8) :: X, sqeps, xbig, y
  INTEGER, EXTERNAL :: INITDS
  REAL(8), EXTERNAL :: DERFC, D1MACH, DCSEVL
  SAVE nterf, xbig, sqeps
  REAL(8), PARAMETER :: erfcs(21) = [ -.49046121234691808039984544033376D-1, &
    -.14226120510371364237824741899631D+0, +.10035582187599795575754676712933D-1, &
    -.57687646997674847650827025509167D-3, +.27419931252196061034422160791471D-4, &
    -.11043175507344507604135381295905D-5, +.38488755420345036949961311498174D-7, &
    -.11808582533875466969631751801581D-8, +.32334215826050909646402930953354D-10, &
    -.79910159470045487581607374708595D-12, +.17990725113961455611967245486634D-13, &
    -.37186354878186926382316828209493D-15, +.71035990037142529711689908394666D-17, &
    -.12612455119155225832495424853333D-18, +.20916406941769294369170500266666D-20, &
    -.32539731029314072982364160000000D-22, +.47668672097976748332373333333333D-24, &
    -.65980120782851343155199999999999D-26, +.86550114699637626197333333333333D-28, &
    -.10788925177498064213333333333333D-29, +.12811883993017002666666666666666D-31 ]
  REAL(8), PARAMETER :: sqrtpi = 1.77245385090551602729816748334115D0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DERF
  IF ( first ) THEN
    nterf = INITDS(erfcs,21,0.1*REAL(D1MACH(3)))
    xbig = SQRT(-LOG(sqrtpi*D1MACH(3)))
    sqeps = SQRT(2.0D0*D1MACH(3))
    first = .FALSE.
  ENDIF
  !
  y = ABS(X)
  IF ( y>1.D0 ) THEN
    !
    ! ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
    !
    IF ( y<=xbig ) THEN
      DERF = SIGN(1.0D0-DERFC(y),X)
    ELSE
      DERF = SIGN(1.0D0,X)
    ENDIF
    RETURN
  ENDIF
  !
  ! ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
  !
  IF ( y<=sqeps ) DERF = 2.0D0*X*X/sqrtpi
  IF ( y>sqeps ) DERF = X*(1.0D0+DCSEVL(2.D0*X*X-1.D0,erfcs,nterf))
  RETURN
END FUNCTION DERF
