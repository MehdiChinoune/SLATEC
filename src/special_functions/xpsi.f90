!** XPSI
REAL(SP) ELEMENTAL FUNCTION XPSI(A,Ipsik,Ipsix)
  !> To compute values of the Psi function for XLEGF.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C7C
  !***
  ! **Type:**      SINGLE PRECISION (XPSI-S, DXPSI-D)
  !***
  ! **Keywords:**  PSI FUNCTION
  !***
  ! **Author:**  Smith, John M., (NBS and George Mason University)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820728  DATE WRITTEN
  !   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
  !   901019  Revisions to prologue.  (DWL and WRB)
  !   901106  Changed all specific intrinsics to generic.  (WRB)
  !           Corrected order of sections in prologue and added TYPE section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER, INTENT(IN) :: Ipsik, Ipsix
  REAL(SP), INTENT(IN) :: A
  INTEGER :: i, k, k1, m, n
  REAL(SP) :: b, c
  !
  !        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
  !        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
  !        NUMBER.
  !
  REAL(SP), PARAMETER :: cnum(12) = [ 1._SP, -1._SP, 1._SP, -1._SP, 1._SP, -691._SP, &
    1._SP, -3617._SP, 43867._SP, -174611._SP, 77683._SP, -236364091._SP ]
  REAL(SP), PARAMETER :: cdenom(12) = [ 12._SP, 120._SP, 252._SP, 240._SP, 132._SP, &
    32760._SP, 12._SP, 8160._SP, 14364._SP, 6600._SP, 276._SP, 65520._SP ]
  !* FIRST EXECUTABLE STATEMENT  XPSI
  n = MAX(0,Ipsix-INT(A))
  b = n + A
  k1 = Ipsik - 1
  !
  !        SERIES EXPANSION FOR A > IPSIX USING IPSIK-1 TERMS.
  !
  c = 0.
  DO i = 1, k1
    k = Ipsik - i
    c = (c+cnum(k)/cdenom(k))/b**2
  END DO
  XPSI = LOG(b) - (c+.5_SP/b)
  IF( n/=0 ) THEN
    b = 0.
    !
    !        RECURRENCE FOR A <= IPSIX.
    !
    DO m = 1, n
      b = b + 1._SP/(n-m+A)
    END DO
    XPSI = XPSI - b
  END IF

END FUNCTION XPSI