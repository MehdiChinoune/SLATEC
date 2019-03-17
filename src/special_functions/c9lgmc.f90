!DECK C9LGMC
COMPLEX FUNCTION C9LGMC(Zin)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  C9LGMC
  !***SUBSIDIARY
  !***PURPOSE  Compute the log gamma correction factor so that
  !            LOG(CGAMMA(Z)) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z
  !            + C9LGMC(Z).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      COMPLEX (R9LGMC-S, D9LGMC-D, C9LGMC-C)
  !***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
  !             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the LOG GAMMA correction term for large ABS(Z) when REAL(Z)
  ! .GE. 0.0 and for large ABS(AIMAG(Y)) when REAL(Z) .LT. 0.0.  We find
  ! C9LGMC so that
  !   LOG(Z) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z)
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  C9LGMC
  REAL bern, bound, cabsz, R1MACH, x, xbig, xmax, y
  INTEGER i, ndx, nterm
  COMPLEX Zin, z, z2inv
  DIMENSION bern(11)
  LOGICAL first
  SAVE bern, nterm, bound, xbig, xmax, first
  DATA bern(1)/.083333333333333333E0/
  DATA bern(2)/ - .0027777777777777778E0/
  DATA bern(3)/.00079365079365079365E0/
  DATA bern(4)/ - .00059523809523809524E0/
  DATA bern(5)/.00084175084175084175E0/
  DATA bern(6)/ - .0019175269175269175E0/
  DATA bern(7)/.0064102564102564103E0/
  DATA bern(8)/ - .029550653594771242E0/
  DATA bern(9)/.17964437236883057E0/
  DATA bern(10)/ - 1.3924322169059011E0/
  DATA bern(11)/13.402864044168392E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  C9LGMC
  IF ( first ) THEN
    nterm = INT( -0.30*LOG(R1MACH(3)) )
    bound = 0.1170*nterm*(0.1*R1MACH(3))**(-1./(2*nterm-1))
    xbig = 1.0/SQRT(R1MACH(3))
    xmax = EXP(MIN(LOG(R1MACH(2)/12.0),-LOG(12.*R1MACH(1))))
  ENDIF
  first = .FALSE.
  !
  z = Zin
  x = REAL(z)
  y = AIMAG(z)
  cabsz = ABS(z)
  !
  IF ( x<0.0.AND.ABS(y)<bound ) CALL XERMSG('SLATEC','C9LGMC',&
    'NOT VALID FOR NEGATIVE REAL(Z) AND SMALL ABS(AIMAG(Z))',2,2)
  IF ( cabsz<bound ) CALL XERMSG('SLATEC','C9LGMC',&
    'NOT VALID FOR SMALL ABS(Z)',3,2)
  !
  IF ( cabsz>=xmax ) THEN
    !
    C9LGMC = (0.0,0.0)
    CALL XERMSG('SLATEC','C9LGMC','Z SO BIG C9LGMC UNDERFLOWS',1,1)
    RETURN
  ENDIF
  !
  IF ( cabsz>=xbig ) C9LGMC = 1.0/(12.0*z)
  IF ( cabsz>=xbig ) RETURN
  !
  z2inv = 1.0/z**2
  C9LGMC = (0.0,0.0)
  DO i = 1, nterm
    ndx = nterm + 1 - i
    C9LGMC = bern(ndx) + C9LGMC*z2inv
  ENDDO
  !
  C9LGMC = C9LGMC/z
  RETURN
END FUNCTION C9LGMC
