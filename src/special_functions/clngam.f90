!** CLNGAM
COMPLEX FUNCTION CLNGAM(Zin)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithm of the absolute value of the Gamma
  !            function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      COMPLEX (ALNGAM-S, DLNGAM-D, CLNGAM-C)
  !***
  ! **Keywords:**  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CLNGAM computes the natural log of the complex valued gamma function
  ! at ZIN, where ZIN is a complex number.  This is a preliminary version,
  ! which is not accurate.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  C9LGMC, CARG, CLNREL, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  
  REAL argsum, bound, cabsz, CARG, dxrel, pi, R1MACH, sq2pil, x, y
  INTEGER i, n
  COMPLEX Zin, z, corr, CLNREL, C9LGMC
  LOGICAL first
  SAVE pi, sq2pil, bound, dxrel, first
  DATA pi/3.14159265358979324E0/
  DATA sq2pil/0.91893853320467274E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  CLNGAM
  IF ( first ) THEN
    n = INT( -0.30*LOG(R1MACH(3)) )
    ! BOUND = N*(0.1*EPS)**(-1/(2*N-1))/(PI*EXP(1))
    bound = 0.1171*n*(0.1*R1MACH(3))**(-1./(2*n-1))
    dxrel = SQRT(R1MACH(4))
  ENDIF
  first = .FALSE.
  !
  z = Zin
  x = REAL(Zin)
  y = AIMAG(Zin)
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
        IF ( x<(-0.5).AND.ABS(y)<=dxrel ) THEN
          IF ( ABS((z-AINT(x-0.5))/x)<dxrel ) CALL XERMSG('SLATEC','CLNGAM',&
            'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER',1,1)
        ENDIF
        !
        n = INT( SQRT(bound**2-y**2) - x ) + 1
        argsum = 0.0
        corr = (1.0,0.0)
        DO i = 1, n
          argsum = argsum + CARG(z)
          corr = z*corr
          z = 1.0 + z
        ENDDO
        !
        IF ( REAL(corr)==0.0.AND.AIMAG(corr)==0.0 )&
          CALL XERMSG('SLATEC','CLNGAM','Z IS A NEGATIVE INTEGER',3,2)
        corr = -CMPLX(LOG(ABS(corr)),argsum)
      ELSE
        !
        ! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
        ! ABS(AIMAG(Y)) SMALL.
        !
        IF ( y>0.0 ) z = CONJG(z)
        corr = EXP(-CMPLX(0.0,2.0*pi)*z)
        IF ( REAL(corr)==1.0.AND.AIMAG(corr)==0.0 )&
          CALL XERMSG('SLATEC','CLNGAM','Z IS A NEGATIVE INTEGER',3,2)
        !
        CLNGAM = sq2pil + 1.0 - CMPLX(0.0,pi)*(z-0.5) - CLNREL(-corr)&
          + (z-0.5)*LOG(1.0-z) - z - C9LGMC(1.0-z)
        IF ( y>0.0 ) CLNGAM = CONJG(CLNGAM)
        RETURN
      ENDIF
    ENDIF
  ENDIF
  !
  ! USE STIRLING-S APPROXIMATION FOR LARGE Z.
  !
  CLNGAM = sq2pil + (z-0.5)*LOG(z) - z + corr + C9LGMC(z)
  !
END FUNCTION CLNGAM
