!DECK CPSI
COMPLEX FUNCTION CPSI(Zin)
  IMPLICIT NONE
  REAL bern, bound, cabsz, dxrel, pi, R1MACH, rbig, rmin, x, y
  INTEGER i, n, ndx, nterm
  !***BEGIN PROLOGUE  CPSI
  !***PURPOSE  Compute the Psi (or Digamma) function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7C
  !***TYPE      COMPLEX (PSI-S, DPSI-D, CPSI-C)
  !***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! PSI(X) calculates the psi (or digamma) function of X.  PSI(X)
  ! is the logarithmic derivative of the gamma function of X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CCOT, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  CPSI
  COMPLEX Zin, z, z2inv, corr, CCOT
  DIMENSION bern(13)
  LOGICAL first
  EXTERNAL CCOT
  SAVE bern, pi, nterm, bound, dxrel, rmin, rbig, first
  DATA bern(1)/.83333333333333333E-1/
  DATA bern(2)/ - .83333333333333333E-2/
  DATA bern(3)/.39682539682539683E-2/
  DATA bern(4)/ - .41666666666666667E-2/
  DATA bern(5)/.75757575757575758E-2/
  DATA bern(6)/ - .21092796092796093E-1/
  DATA bern(7)/.83333333333333333E-1/
  DATA bern(8)/ - .44325980392156863E0/
  DATA bern(9)/.30539543302701197E1/
  DATA bern(10)/ - .26456212121212121E2/
  DATA bern(11)/.28146014492753623E3/
  DATA bern(12)/ - .34548853937728938E4/
  DATA bern(13)/.54827583333333333E5/
  DATA pi/3.141592653589793E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  CPSI
  IF ( first ) THEN
    nterm = -0.30*LOG(R1MACH(3))
    ! MAYBE BOUND = N*(0.1*EPS)**(-1/(2*N-1)) / (PI*EXP(1))
    bound = 0.1171*nterm*(0.1*R1MACH(3))**(-1.0/(2*nterm-1))
    dxrel = SQRT(R1MACH(4))
    rmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+0.011)
    rbig = 1.0/R1MACH(3)
  ENDIF
  first = .FALSE.
  !
  z = Zin
  x = REAL(z)
  y = AIMAG(z)
  IF ( y<0.0 ) z = CONJG(z)
  !
  corr = (0.0,0.0)
  cabsz = ABS(z)
  IF ( x<0.0.OR.cabsz<=bound ) THEN
    IF ( x>=0.0.OR.ABS(y)<=bound ) THEN
      !
      IF ( cabsz<bound ) THEN
        !
        ! USE THE RECURSION RELATION FOR ABS(Z) SMALL.
        !
        IF ( cabsz<rmin ) CALL XERMSG('SLATEC','CPSI',&
          'CPSI CALLED WITH Z SO NEAR 0 THAT CPSI OVERFLOWS',2,2)
        !
        IF ( x<(-0.5).AND.ABS(y)<=dxrel ) THEN
          IF ( ABS((z-AINT(x-0.5))/x)<dxrel ) CALL XERMSG('SLATEC','CPSI',&
            'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER'&
            ,1,1)
          IF ( y==0.0.AND.x==AINT(x) )&
            CALL XERMSG('SLATEC','CPSI','Z IS A NEGATIVE INTEGER',3,2)
        ENDIF
        !
        n = SQRT(bound**2-y**2) - x + 1.0
        DO i = 1, n
          corr = corr - 1.0/z
          z = z + 1.0
        ENDDO
      ELSE
        !
        ! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
        ! ABS(AIMAG(Y)) SMALL.
        !
        corr = -pi*CCOT(pi*z)
        z = 1.0 - z
      ENDIF
    ENDIF
  ENDIF
  !
  ! NOW EVALUATE THE ASYMPTOTIC SERIES FOR SUITABLY LARGE Z.
  !
  IF ( cabsz>rbig ) CPSI = LOG(z) + corr
  IF ( cabsz<=rbig ) THEN
    !
    CPSI = (0.0,0.0)
    z2inv = 1.0/z**2
    DO i = 1, nterm
      ndx = nterm + 1 - i
      CPSI = bern(ndx) + z2inv*CPSI
    ENDDO
    CPSI = LOG(z) - 0.5/z - CPSI*z2inv + corr
  ENDIF
  !
  IF ( y<0.0 ) CPSI = CONJG(CPSI)
  !
END FUNCTION CPSI
