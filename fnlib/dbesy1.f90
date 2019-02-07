!*==DBESY1.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBESY1
DOUBLE PRECISION FUNCTION DBESY1(X)
  IMPLICIT NONE
  !*--DBESY15
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS , nty1
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBESY1
  !***PURPOSE  Compute the Bessel function of the second kind of order
  !            one.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10A1
  !***TYPE      DOUBLE PRECISION (BESY1-S, DBESY1-D)
  !***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DBESY1(X) calculates the double precision Bessel function of the
  ! second kind of order for double precision argument X.
  !
  ! Series for BY1        on the interval  0.          to  1.60000E+01
  !                                        with weighted error   8.65E-33
  !                                         log weighted error  32.06
  !                               significant figures required  32.17
  !                                    decimal places required  32.71
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9B1MP, DBESJ1, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBESY1
  DOUBLE PRECISION X , by1cs(20) , ampl , theta , twodpi , xmin , xsml , y , &
    D1MACH , DCSEVL , DBESJ1
  LOGICAL first
  SAVE by1cs , twodpi , nty1 , xmin , xsml , first
  DATA by1cs(1)/ + .320804710061190862932352018628015D-1/
  DATA by1cs(2)/ + .126270789743350044953431725999727D+1/
  DATA by1cs(3)/ + .649996189992317500097490637314144D-2/
  DATA by1cs(4)/ - .893616452886050411653144160009712D-1/
  DATA by1cs(5)/ + .132508812217570954512375510370043D-1/
  DATA by1cs(6)/ - .897905911964835237753039508298105D-3/
  DATA by1cs(7)/ + .364736148795830678242287368165349D-4/
  DATA by1cs(8)/ - .100137438166600055549075523845295D-5/
  DATA by1cs(9)/ + .199453965739017397031159372421243D-7/
  DATA by1cs(10)/ - .302306560180338167284799332520743D-9/
  DATA by1cs(11)/ + .360987815694781196116252914242474D-11/
  DATA by1cs(12)/ - .348748829728758242414552947409066D-13/
  DATA by1cs(13)/ + .278387897155917665813507698517333D-15/
  DATA by1cs(14)/ - .186787096861948768766825352533333D-17/
  DATA by1cs(15)/ + .106853153391168259757070336000000D-19/
  DATA by1cs(16)/ - .527472195668448228943872000000000D-22/
  DATA by1cs(17)/ + .227019940315566414370133333333333D-24/
  DATA by1cs(18)/ - .859539035394523108693333333333333D-27/
  DATA by1cs(19)/ + .288540437983379456000000000000000D-29/
  DATA by1cs(20)/ - .864754113893717333333333333333333D-32/
  DATA twodpi/0.636619772367581343075535053490057D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBESY1
  IF ( first ) THEN
    nty1 = INITDS(by1cs,20,0.1*REAL(D1MACH(3)))
    !
    xmin = 1.571D0*EXP(MAX(LOG(D1MACH(1)),-LOG(D1MACH(2)))+0.01D0)
    xsml = SQRT(4.0D0*D1MACH(3))
  ENDIF
  first = .FALSE.
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DBESY1','X IS ZERO OR NEGATIVE',1,2)
  IF ( X>4.0D0 ) THEN
    !
    CALL D9B1MP(X,ampl,theta)
    DBESY1 = ampl*SIN(theta)
    GOTO 99999
  ENDIF
  !
  IF ( X<xmin ) CALL XERMSG('SLATEC','DBESY1','X SO SMALL Y1 OVERFLOWS',3,2)
  y = 0.D0
  IF ( X>xsml ) y = X*X
  DBESY1 = twodpi*LOG(0.5D0*X)*DBESJ1(X)&
    + (0.5D0+DCSEVL(.125D0*y-1.D0,by1cs,nty1))/X
  RETURN
  !
  99999 END FUNCTION DBESY1
