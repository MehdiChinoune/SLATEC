!*==ASINH.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK ASINH
FUNCTION ASINH(X)
  IMPLICIT NONE
  !*--ASINH5
  !*** Start of declarations inserted by SPAG
  REAL aln2, ASINH, asnhcs, CSEVL, R1MACH, sqeps, X, xmax, y
  INTEGER INITS, nterms
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  ASINH
  !***PURPOSE  Compute the arc hyperbolic sine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      SINGLE PRECISION (ASINH-S, DASINH-D, CASINH-C)
  !***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC SINE
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ASINH(X) computes the arc hyperbolic sine of X.
  !
  ! Series for ASNH       on the interval  0.          to  1.00000D+00
  !                                        with weighted error   2.19E-17
  !                                         log weighted error  16.66
  !                               significant figures required  15.60
  !                                    decimal places required  17.31
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  ASINH
  DIMENSION asnhcs(20)
  LOGICAL first
  SAVE aln2, asnhcs, nterms, xmax, sqeps, first
  DATA aln2/0.69314718055994530942E0/
  DATA asnhcs(1)/ - .12820039911738186E0/
  DATA asnhcs(2)/ - .058811761189951768E0/
  DATA asnhcs(3)/.004727465432212481E0/
  DATA asnhcs(4)/ - .000493836316265361E0/
  DATA asnhcs(5)/.000058506207058557E0/
  DATA asnhcs(6)/ - .000007466998328931E0/
  DATA asnhcs(7)/.000001001169358355E0/
  DATA asnhcs(8)/ - .000000139035438587E0/
  DATA asnhcs(9)/.000000019823169483E0/
  DATA asnhcs(10)/ - .000000002884746841E0/
  DATA asnhcs(11)/.000000000426729654E0/
  DATA asnhcs(12)/ - .000000000063976084E0/
  DATA asnhcs(13)/.000000000009699168E0/
  DATA asnhcs(14)/ - .000000000001484427E0/
  DATA asnhcs(15)/.000000000000229037E0/
  DATA asnhcs(16)/ - .000000000000035588E0/
  DATA asnhcs(17)/.000000000000005563E0/
  DATA asnhcs(18)/ - .000000000000000874E0/
  DATA asnhcs(19)/.000000000000000138E0/
  DATA asnhcs(20)/ - .000000000000000021E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  ASINH
  IF ( first ) THEN
    nterms = INITS(asnhcs,20,0.1*R1MACH(3))
    sqeps = SQRT(R1MACH(3))
    xmax = 1.0/sqeps
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>1.0 ) THEN
    !
    IF ( y<xmax ) ASINH = LOG(y+SQRT(y**2+1.))
    IF ( y>=xmax ) ASINH = aln2 + LOG(y)
    ASINH = SIGN(ASINH,X)
    GOTO 99999
  ENDIF
  !
  ASINH = X
  IF ( y>sqeps ) ASINH = X*(1.0+CSEVL(2.*X*X-1.,asnhcs,nterms))
  RETURN
  !
  99999 CONTINUE
  END FUNCTION ASINH
