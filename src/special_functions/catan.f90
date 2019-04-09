!** CATAN
COMPLEX FUNCTION CATAN(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex arc tangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      COMPLEX (CATAN-C)
  !***
  ! **Keywords:**  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CATAN(Z) calculates the complex trigonometric arc tangent of Z.
  ! The result is in units of radians, and the real part is in the first
  ! or fourth quadrant.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  INTEGER i
  REAL r, R1MACH, r2, twoi, x, xans, y, yans
  COMPLEX Z, z2
  INTEGER, SAVE :: nterms
  REAL, SAVE :: sqeps, rmin, rmax
  REAL, PARAMETER :: pi2 = 1.57079632679489661923E0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  CATAN
  IF ( first ) THEN
    ! NTERMS = LOG(EPS)/LOG(RBND) WHERE RBND = 0.1
    nterms = INT( -0.4343*LOG(R1MACH(3)) ) + 1
    sqeps = SQRT(R1MACH(4))
    rmin = SQRT(3.0*R1MACH(3))
    rmax = 1.0/R1MACH(3)
    first = .FALSE.
  END IF
  !
  r = ABS(Z)
  IF ( r<=0.1 ) THEN
    !
    CATAN = Z
    IF ( r<rmin ) RETURN
    !
    CATAN = (0.0,0.0)
    z2 = Z*Z
    DO i = 1, nterms
      twoi = 2*(nterms-i) + 1
      CATAN = 1.0/twoi - z2*CATAN
    END DO
    CATAN = Z*CATAN
    RETURN
    !
  ELSEIF ( r>rmax ) THEN
    !
    CATAN = CMPLX(pi2,0.)
    IF ( REAL(Z)<0.0 ) CATAN = CMPLX(-pi2,0.0)
    RETURN
  ELSE
    x = REAL(Z)
    IF( ABS(x) == 0. ) x = 0.
    y = AIMAG(Z)
    r2 = r*r
    IF ( r2==1.0.AND.x==0.0 ) CALL XERMSG('SLATEC','CATAN','Z IS +I OR -I',2,2)
    IF ( ABS(r2-1.0)<=sqeps ) THEN
      IF ( ABS(CMPLX(1.0,0.0)+Z*Z)<sqeps ) CALL XERMSG('SLATEC','CATAN',&
        'ANSWER LT HALF PRECISION, Z**2 CLOSE TO -1',1,1)
    END IF
  END IF
  !
  xans = 0.5*ATAN2(2.0*x,1.0-r2)
  yans = 0.25*LOG((r2+2.0*y+1.0)/(r2-2.0*y+1.0))
  CATAN = CMPLX(xans,yans)
  RETURN
END FUNCTION CATAN
