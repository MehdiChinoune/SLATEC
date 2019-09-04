!** CLNGAM
COMPLEX(SP) ELEMENTAL FUNCTION CLNGAM(Zin)
  !> Compute the logarithm of the absolute value of the Gamma function.
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
  USE service, ONLY : eps_2_sp
  !
  COMPLEX(SP), INTENT(IN) :: Zin
  !
  INTEGER :: i, n
  REAL(SP) :: argsum, cabsz, x, y
  COMPLEX(SP) :: z, corr
  INTEGER, PARAMETER :: np = INT( -0.30_SP*LOG(eps_2_sp) )
  ! BOUND = N*(0.1*EPS)**(-1/(2*N-1))/(PI*EXP(1))
  REAL(SP), PARAMETER :: bound = 0.1171_SP*np*(0.1_SP*eps_2_sp)**(-1._SP/(2*np-1))
  REAL(SP), PARAMETER :: pi = 3.14159265358979324_SP
  REAL(SP), PARAMETER :: sq2pil = 0.91893853320467274_SP
  !* FIRST EXECUTABLE STATEMENT  CLNGAM
  !
  z = Zin
  x = REAL(Zin)
  y = AIMAG(Zin)
  !
  corr = (0._SP,0._SP)
  cabsz = ABS(z)
  IF( x<0._SP .OR. cabsz<=bound ) THEN
    IF( x>=0._SP .OR. ABS(y)<=bound ) THEN
      !
      IF( cabsz<bound ) THEN
        !
        ! USE THE RECURSION RELATION FOR ABS(Z) SMALL.
        !
        ! IF( x<(-0.5_SP) .AND. ABS(y)<=dxrel ) THEN
        !   IF( ABS((z-AINT(x-0.5_SP))/x)<dxrel ) CALL XERMSG('CLNGAM',&
        !     'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER',1,1)
        ! END IF
        !
        n = INT( SQRT(bound**2-y**2) - x ) + 1
        argsum = 0._SP
        corr = (1._SP,0._SP)
        DO i = 1, n
          argsum = argsum + CARG(z)
          corr = z*corr
          z = 1._SP + z
        END DO
        !
        IF( REAL(corr)==0._SP .AND. AIMAG(corr)==0._SP ) THEN
          ERROR STOP 'CLNGAM : Z IS A NEGATIVE INTEGER'
        END IF
        corr = -CMPLX(LOG(ABS(corr)),argsum,SP)
      ELSE
        !
        ! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
        ! ABS(AIMAG(Y)) SMALL.
        !
        IF( y>0._SP ) z = CONJG(z)
        corr = EXP(-CMPLX(0._SP,2._SP*pi,SP)*z)
        IF( REAL(corr)==1._SP .AND. AIMAG(corr)==0._SP ) THEN
          ERROR STOP 'CLNGAM : Z IS A NEGATIVE INTEGER'
        END IF
        !
        CLNGAM = sq2pil + 1._SP - CMPLX(0._SP,pi,SP)*(z-0.5_SP) - CLNREL(-corr)&
          + (z-0.5_SP)*LOG(1._SP-z) - z - C9LGMC(1._SP-z)
        IF( y>0._SP ) CLNGAM = CONJG(CLNGAM)
        RETURN
      END IF
    END IF
  END IF
  !
  ! USE STIRLING-S APPROXIMATION FOR LARGE Z.
  !
  CLNGAM = sq2pil + (z-0.5_SP)*LOG(z) - z + corr + C9LGMC(z)
  !
END FUNCTION CLNGAM