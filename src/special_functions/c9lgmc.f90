!** C9LGMC
COMPLEX FUNCTION C9LGMC(Zin)
  !>
  !***
  !  Compute the log gamma correction factor so that
  !            LOG(CGAMMA(Z)) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z
  !            + C9LGMC(Z).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      COMPLEX (R9LGMC-S, D9LGMC-D, C9LGMC-C)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
  !             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the LOG GAMMA correction term for large ABS(Z) when REAL(Z)
  ! .GE. 0.0 and for large ABS(AIMAG(Y)) when REAL(Z) .LT. 0.0.  We find
  ! C9LGMC so that
  !   LOG(Z) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z)
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)

  REAL cabsz, x, y
  INTEGER i, ndx
  COMPLEX Zin, z, z2inv
  INTEGER, SAVE :: nterm
  REAL, SAVE :: bound, xbig, xmax
  REAL, PARAMETER :: bern(11) = [ .083333333333333333E0,-.0027777777777777778E0, &
    .00079365079365079365E0, -.00059523809523809524E0, .00084175084175084175E0, &
    -.0019175269175269175E0,  .0064102564102564103E0, -.029550653594771242E0, &
    .17964437236883057E0, -1.3924322169059011E0,   13.402864044168392E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  C9LGMC
  IF ( first ) THEN
    nterm = INT( -0.30*LOG(R1MACH(3)) )
    bound = 0.1170*nterm*(0.1*R1MACH(3))**(-1./(2*nterm-1))
    xbig = 1.0/SQRT(R1MACH(3))
    xmax = EXP(MIN(LOG(R1MACH(2)/12.0),-LOG(12.*R1MACH(1))))
    first = .FALSE.
  END IF
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
  END IF
  !
  IF ( cabsz>=xbig ) C9LGMC = 1.0/(12.0*z)
  IF ( cabsz>=xbig ) RETURN
  !
  z2inv = 1.0/z**2
  C9LGMC = (0.0,0.0)
  DO i = 1, nterm
    ndx = nterm + 1 - i
    C9LGMC = bern(ndx) + C9LGMC*z2inv
  END DO
  !
  C9LGMC = C9LGMC/z
  RETURN
END FUNCTION C9LGMC
