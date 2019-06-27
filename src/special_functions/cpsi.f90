!** CPSI
COMPLEX(SP) ELEMENTAL FUNCTION CPSI(Zin)
  !> Compute the Psi (or Digamma) function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7C
  !***
  ! **Type:**      COMPLEX (PSI-S, DPSI-D, CPSI-C)
  !***
  ! **Keywords:**  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! PSI(X) calculates the psi (or digamma) function of X.  PSI(X)
  ! is the logarithmic derivative of the gamma function of X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CCOT, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780501  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  USE service, ONLY : R1MACH
  COMPLEX(SP), INTENT(IN) :: Zin
  INTEGER :: i, n, ndx
  REAL(SP) :: cabsz, x, y
  COMPLEX(SP) :: z, z2inv, corr
  INTEGER, PARAMETER :: nterm = INT( -0.30*LOG(R1MACH(3)) )
  ! MAYBE BOUND = N*(0.1*EPS)**(-1/(2*N-1)) / (PI*EXP(1))
  REAL(SP), PARAMETER ::  bound = 0.1171_SP*nterm*(0.1_SP*R1MACH(3))**(-1._SP/(2*nterm-1)), &
    dxrel = SQRT(R1MACH(4)), rmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+0.011_SP), &
    rbig = 1._SP/R1MACH(3)
  REAL(SP), PARAMETER :: bern(13) = [ .83333333333333333E-1_SP, -.83333333333333333E-2_SP, &
    .39682539682539683E-2_SP, -.41666666666666667E-2_SP, .75757575757575758E-2_SP, &
    -.21092796092796093E-1_SP, .83333333333333333E-1_SP, -.44325980392156863E0_SP, &
    .30539543302701197E1_SP,  -.26456212121212121E2_SP,  .28146014492753623E3_SP, &
    -.34548853937728938E4_SP,  .54827583333333333_SP ]
  REAL(SP), PARAMETER :: pi = 3.141592653589793_SP
  !* FIRST EXECUTABLE STATEMENT  CPSI
  !
  z = Zin
  x = REAL(z)
  y = AIMAG(z)
  IF( y<0._SP ) z = CONJG(z)
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
        IF( cabsz<rmin ) ERROR STOP 'CPSI : CPSI CALLED WITH Z SO NEAR 0 THAT CPSI OVERFLOWS'
        IF( y==0._SP .AND. x==AINT(x) ) ERROR STOP 'CPSI : Z IS A NEGATIVE INTEGER'
        !
        ! IF( x<(-0.5_SP) .AND. ABS(y)<=dxrel .AND. ABS((z-AINT(x-0.5_SP))/x)<dxrel ) THEN
          ! CALL XERMSG('CPSI : ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER'
        ! END IF
        !
        n = INT( SQRT(bound**2-y**2) - x ) + 1
        DO i = 1, n
          corr = corr - 1._SP/z
          z = z + 1._SP
        END DO
      ELSE
        !
        ! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
        ! ABS(AIMAG(Y)) SMALL.
        !
        corr = -pi*CCOT(pi*z)
        z = 1._SP - z
      END IF
    END IF
  END IF
  !
  ! NOW EVALUATE THE ASYMPTOTIC SERIES FOR SUITABLY LARGE Z.
  !
  IF( cabsz>rbig ) CPSI = LOG(z) + corr
  IF( cabsz<=rbig ) THEN
    !
    CPSI = (0._SP,0._SP)
    z2inv = 1._SP/z**2
    DO i = 1, nterm
      ndx = nterm + 1 - i
      CPSI = bern(ndx) + z2inv*CPSI
    END DO
    CPSI = LOG(z) - 0.5_SP/z - CPSI*z2inv + corr
  END IF
  !
  IF( y<0._SP ) CPSI = CONJG(CPSI)
  !
END FUNCTION CPSI