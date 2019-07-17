!** C9LGMC
COMPLEX(SP) ELEMENTAL FUNCTION C9LGMC(Zin)
  !> Compute the log gamma correction factor so that
  !    LOG(CGAMMA(Z)) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z).
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
  ! Compute the LOG GAMMA correction term for large ABS(Z) when REAL(Z) >= 0.0
  ! and for large ABS(AIMAG(Y)) when REAL(Z) < 0.0.  We find C9LGMC so that
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
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : R1MACH
  COMPLEX(SP), INTENT(IN) :: Zin
  INTEGER :: i, ndx
  REAL(SP) :: cabsz, x, y
  COMPLEX(SP) :: z, z2inv
  INTEGER, PARAMETER :: nterm = INT( -0.30*LOG(R1MACH(3)) )
  REAL(SP), PARAMETER :: bound = 0.117_SP*nterm*(0.1_SP*R1MACH(3))**(-1._SP/(2*nterm-1)), &
    xbig = 1._SP/SQRT(R1MACH(3)), &
    xmax = EXP(MIN(LOG(R1MACH(2)/12._SP),-LOG(12._SP*R1MACH(1))))
  REAL(SP), PARAMETER :: bern(11) = [ .083333333333333333_SP,-.0027777777777777778_SP, &
    .00079365079365079365_SP, -.00059523809523809524_SP, .00084175084175084175_SP, &
    -.0019175269175269175_SP,  .0064102564102564103_SP, -.029550653594771242_SP, &
    .17964437236883057_SP, -1.3924322169059011_SP,   13.402864044168392_SP ]
  !* FIRST EXECUTABLE STATEMENT  C9LGMC
  !
  z = Zin
  x = REAL(z)
  y = AIMAG(z)
  cabsz = ABS(z)
  !
  IF( x<0._SP .AND. ABS(y)<bound ) THEN
    ERROR STOP 'C9LGMC : &NOT VALID FOR NEGATIVE REAL(Z) AND SMALL ABS(AIMAG(Z))'
  ELSEIF( cabsz<bound ) THEN
    ERROR STOP 'C9LGMC : NOT VALID FOR SMALL ABS(Z)'
  ELSEIF( cabsz>=xmax ) THEN
    C9LGMC = (0._SP,0._SP)
    ! CALL XERMSG('C9LGMC : Z SO BIG C9LGMC UNDERFLOWS',1,1)
  ELSEIF( cabsz>=xbig ) THEN
    C9LGMC = 1._SP/(12._SP*z)
  ELSE
    z2inv = 1._SP/z**2
    C9LGMC = (0._SP,0._SP)
    DO i = 1, nterm
      ndx = nterm + 1 - i
      C9LGMC = bern(ndx) + C9LGMC*z2inv
    END DO
    C9LGMC = C9LGMC/z
  END IF

  RETURN
END FUNCTION C9LGMC