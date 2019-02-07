!*==DERF.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DERF
REAL(8) FUNCTION DERF(X)
  IMPLICIT NONE
  !*--DERF5
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , nterf
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DERF
  !***PURPOSE  Compute the error function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C8A, L5A1E
  !***TYPE      DOUBLE PRECISION (ERF-S, DERF-D)
  !***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, DERFC, INITDS
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !   920618  Removed space from variable name.  (RWC, WRB)
  !***END PROLOGUE  DERF
  REAL(8) :: X , erfcs(21) , sqeps , sqrtpi , xbig , y , D1MACH , &
    DCSEVL , DERFC
  LOGICAL first
  EXTERNAL DERFC
  SAVE erfcs , sqrtpi , nterf , xbig , sqeps , first
  DATA erfcs(1)/ - .49046121234691808039984544033376D-1/
  DATA erfcs(2)/ - .14226120510371364237824741899631D+0/
  DATA erfcs(3)/ + .10035582187599795575754676712933D-1/
  DATA erfcs(4)/ - .57687646997674847650827025509167D-3/
  DATA erfcs(5)/ + .27419931252196061034422160791471D-4/
  DATA erfcs(6)/ - .11043175507344507604135381295905D-5/
  DATA erfcs(7)/ + .38488755420345036949961311498174D-7/
  DATA erfcs(8)/ - .11808582533875466969631751801581D-8/
  DATA erfcs(9)/ + .32334215826050909646402930953354D-10/
  DATA erfcs(10)/ - .79910159470045487581607374708595D-12/
  DATA erfcs(11)/ + .17990725113961455611967245486634D-13/
  DATA erfcs(12)/ - .37186354878186926382316828209493D-15/
  DATA erfcs(13)/ + .71035990037142529711689908394666D-17/
  DATA erfcs(14)/ - .12612455119155225832495424853333D-18/
  DATA erfcs(15)/ + .20916406941769294369170500266666D-20/
  DATA erfcs(16)/ - .32539731029314072982364160000000D-22/
  DATA erfcs(17)/ + .47668672097976748332373333333333D-24/
  DATA erfcs(18)/ - .65980120782851343155199999999999D-26/
  DATA erfcs(19)/ + .86550114699637626197333333333333D-28/
  DATA erfcs(20)/ - .10788925177498064213333333333333D-29/
  DATA erfcs(21)/ + .12811883993017002666666666666666D-31/
  DATA sqrtpi/1.77245385090551602729816748334115D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DERF
  IF ( first ) THEN
    nterf = INITDS(erfcs,21,0.1*REAL(D1MACH(3)))
    xbig = SQRT(-LOG(sqrtpi*D1MACH(3)))
    sqeps = SQRT(2.0D0*D1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1.D0 ) THEN
    !
    ! ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
    !
    IF ( y<=xbig ) DERF = SIGN(1.0D0-DERFC(y),X)
    IF ( y>xbig ) DERF = SIGN(1.0D0,X)
    GOTO 99999
  ENDIF
  !
  ! ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
  !
  IF ( y<=sqeps ) DERF = 2.0D0*X*X/sqrtpi
  IF ( y>sqeps ) DERF = X*(1.0D0+DCSEVL(2.D0*X*X-1.D0,erfcs,nterf))
  RETURN
  !
  99999 END FUNCTION DERF
