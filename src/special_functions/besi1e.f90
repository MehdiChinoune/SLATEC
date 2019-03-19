!** BESI1E
REAL FUNCTION BESI1E(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the exponentially scaled modified (hyperbolic)
  !            Bessel function of the first kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESI1E-S, DBSI1E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
  !             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
  !             ORDER ONE, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESI1E(X) calculates the exponentially scaled modified (hyperbolic)
  ! Bessel function of the first kind of order one for real argument X;
  ! i.e., EXP(-ABS(X))*I1(X).
  !
  ! Series for BI1        on the interval  0.          to  9.00000D+00
  !                                        with weighted error   2.40E-17
  !                                         log weighted error  16.62
  !                               significant figures required  16.23
  !                                    decimal places required  17.14
  !
  ! Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
  !                                        with weighted error   6.98E-17
  !                                         log weighted error  16.16
  !                               significant figures required  14.53
  !                                    decimal places required  16.82
  !
  ! Series for AI12       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   3.55E-17
  !                                         log weighted error  16.45
  !                               significant figures required  14.69
  !                                    decimal places required  17.12
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890210  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  
  REAL ai12cs, ai1cs, bi1cs, CSEVL, R1MACH, X, xmin, xsml, y
  INTEGER INITS, ntai1, ntai12, nti1
  DIMENSION bi1cs(11), ai1cs(21), ai12cs(22)
  LOGICAL first
  SAVE bi1cs, ai1cs, ai12cs, nti1, ntai1, ntai12, xmin, xsml, first
  DATA bi1cs(1)/ - .001971713261099859E0/
  DATA bi1cs(2)/.40734887667546481E0/
  DATA bi1cs(3)/.034838994299959456E0/
  DATA bi1cs(4)/.001545394556300123E0/
  DATA bi1cs(5)/.000041888521098377E0/
  DATA bi1cs(6)/.000000764902676483E0/
  DATA bi1cs(7)/.000000010042493924E0/
  DATA bi1cs(8)/.000000000099322077E0/
  DATA bi1cs(9)/.000000000000766380E0/
  DATA bi1cs(10)/.000000000000004741E0/
  DATA bi1cs(11)/.000000000000000024E0/
  DATA ai1cs(1)/ - .02846744181881479E0/
  DATA ai1cs(2)/ - .01922953231443221E0/
  DATA ai1cs(3)/ - .00061151858579437E0/
  DATA ai1cs(4)/ - .00002069971253350E0/
  DATA ai1cs(5)/.00000858561914581E0/
  DATA ai1cs(6)/.00000104949824671E0/
  DATA ai1cs(7)/ - .00000029183389184E0/
  DATA ai1cs(8)/ - .00000001559378146E0/
  DATA ai1cs(9)/.00000001318012367E0/
  DATA ai1cs(10)/ - .00000000144842341E0/
  DATA ai1cs(11)/ - .00000000029085122E0/
  DATA ai1cs(12)/.00000000012663889E0/
  DATA ai1cs(13)/ - .00000000001664947E0/
  DATA ai1cs(14)/ - .00000000000166665E0/
  DATA ai1cs(15)/.00000000000124260E0/
  DATA ai1cs(16)/ - .00000000000027315E0/
  DATA ai1cs(17)/.00000000000002023E0/
  DATA ai1cs(18)/.00000000000000730E0/
  DATA ai1cs(19)/ - .00000000000000333E0/
  DATA ai1cs(20)/.00000000000000071E0/
  DATA ai1cs(21)/ - .00000000000000006E0/
  DATA ai12cs(1)/.02857623501828014E0/
  DATA ai12cs(2)/ - .00976109749136147E0/
  DATA ai12cs(3)/ - .00011058893876263E0/
  DATA ai12cs(4)/ - .00000388256480887E0/
  DATA ai12cs(5)/ - .00000025122362377E0/
  DATA ai12cs(6)/ - .00000002631468847E0/
  DATA ai12cs(7)/ - .00000000383538039E0/
  DATA ai12cs(8)/ - .00000000055897433E0/
  DATA ai12cs(9)/ - .00000000001897495E0/
  DATA ai12cs(10)/.00000000003252602E0/
  DATA ai12cs(11)/.00000000001412580E0/
  DATA ai12cs(12)/.00000000000203564E0/
  DATA ai12cs(13)/ - .00000000000071985E0/
  DATA ai12cs(14)/ - .00000000000040836E0/
  DATA ai12cs(15)/ - .00000000000002101E0/
  DATA ai12cs(16)/.00000000000004273E0/
  DATA ai12cs(17)/.00000000000001041E0/
  DATA ai12cs(18)/ - .00000000000000382E0/
  DATA ai12cs(19)/ - .00000000000000186E0/
  DATA ai12cs(20)/.00000000000000033E0/
  DATA ai12cs(21)/.00000000000000028E0/
  DATA ai12cs(22)/ - .00000000000000003E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  BESI1E
  IF ( first ) THEN
    nti1 = INITS(bi1cs,11,0.1*R1MACH(3))
    ntai1 = INITS(ai1cs,21,0.1*R1MACH(3))
    ntai12 = INITS(ai12cs,22,0.1*R1MACH(3))
    !
    xmin = 2.0*R1MACH(1)
    xsml = SQRT(4.5*R1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>3.0 ) THEN
    !
    IF ( y<=8. ) BESI1E = (.375+CSEVL((48./y-11.)/5.,ai1cs,ntai1))/SQRT(y)
    IF ( y>8. ) BESI1E = (.375+CSEVL(16./y-1.0,ai12cs,ntai12))/SQRT(y)
    BESI1E = SIGN(BESI1E,X)
    RETURN
  ENDIF
  !
  BESI1E = 0.0
  IF ( y==0.0 ) RETURN
  !
  IF ( y<=xmin ) CALL XERMSG('SLATEC','BESI1E',&
    'ABS(X) SO SMALL I1 UNDERFLOWS',1,1)
  IF ( y>xmin ) BESI1E = 0.5*X
  IF ( y>xsml ) BESI1E = X*(.875+CSEVL(y*y/4.5-1.,bi1cs,nti1))
  BESI1E = EXP(-y)*BESI1E
  RETURN
END FUNCTION BESI1E
