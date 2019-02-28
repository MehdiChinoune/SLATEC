!DECK DBESJ1
REAL(8) FUNCTION DBESJ1(X)
  IMPLICIT NONE
  INTEGER INITDS, ntj1
  !***BEGIN PROLOGUE  DBESJ1
  !***PURPOSE  Compute the Bessel function of the first kind of order one.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10A1
  !***TYPE      DOUBLE PRECISION (BESJ1-S, DBESJ1-D)
  !***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBESJ1(X) calculates the double precision Bessel function of the
  ! first kind of order one for double precision argument X.
  !
  ! Series for BJ1        on the interval  0.          to  1.60000E+01
  !                                        with weighted error   1.16E-33
  !                                         log weighted error  32.93
  !                               significant figures required  32.36
  !                                    decimal places required  33.57
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9B1MP, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   910401  Corrected error in code which caused values to have the
  !           wrong sign for arguments less than 4.0.  (WRB)
  !***END PROLOGUE  DBESJ1
  REAL(8) :: X, bj1cs(19), ampl, theta, xsml, xmin, y, D1MACH, &
    DCSEVL
  LOGICAL first
  SAVE bj1cs, ntj1, xsml, xmin, first
  DATA bj1cs(1)/ - .117261415133327865606240574524003D+0/
  DATA bj1cs(2)/ - .253615218307906395623030884554698D+0/
  DATA bj1cs(3)/ + .501270809844695685053656363203743D-1/
  DATA bj1cs(4)/ - .463151480962508191842619728789772D-2/
  DATA bj1cs(5)/ + .247996229415914024539124064592364D-3/
  DATA bj1cs(6)/ - .867894868627882584521246435176416D-5/
  DATA bj1cs(7)/ + .214293917143793691502766250991292D-6/
  DATA bj1cs(8)/ - .393609307918317979229322764073061D-8/
  DATA bj1cs(9)/ + .559118231794688004018248059864032D-10/
  DATA bj1cs(10)/ - .632761640466139302477695274014880D-12/
  DATA bj1cs(11)/ + .584099161085724700326945563268266D-14/
  DATA bj1cs(12)/ - .448253381870125819039135059199999D-16/
  DATA bj1cs(13)/ + .290538449262502466306018688000000D-18/
  DATA bj1cs(14)/ - .161173219784144165412118186666666D-20/
  DATA bj1cs(15)/ + .773947881939274637298346666666666D-23/
  DATA bj1cs(16)/ - .324869378211199841143466666666666D-25/
  DATA bj1cs(17)/ + .120223767722741022720000000000000D-27/
  DATA bj1cs(18)/ - .395201221265134933333333333333333D-30/
  DATA bj1cs(19)/ + .116167808226645333333333333333333D-32/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBESJ1
  IF ( first ) THEN
    ntj1 = INITDS(bj1cs,19,0.1*REAL(D1MACH(3)))
    !
    xsml = SQRT(8.0D0*D1MACH(3))
    xmin = 2.0D0*D1MACH(1)
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>4.0D0 ) THEN
    !
    CALL D9B1MP(y,ampl,theta)
    DBESJ1 = SIGN(ampl,X)*COS(theta)
    RETURN
  ENDIF
  !
  DBESJ1 = 0.0D0
  IF ( y==0.0D0 ) RETURN
  IF ( y<=xmin ) CALL XERMSG('SLATEC','DBESJ1',&
    'ABS(X) SO SMALL J1 UNDERFLOWS',1,1)
  IF ( y>xmin ) DBESJ1 = 0.5D0*X
  IF ( y>xsml ) DBESJ1 = X*(.25D0+DCSEVL(.125D0*y*y-1.D0,bj1cs,ntj1))
  RETURN
END FUNCTION DBESJ1
