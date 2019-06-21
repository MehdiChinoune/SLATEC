!** DXPSI
REAL(DP) FUNCTION DXPSI(A,Ipsik,Ipsix)
  !> To compute values of the Psi function for DXLEGF.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C7C
  !***
  ! **Type:**      DOUBLE PRECISION (XPSI-S, DXPSI-D)
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
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER :: i, Ipsik, Ipsix, k, k1, m, n
  REAL(DP) :: A, b, c
  !
  !        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
  !        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
  !        NUMBER.
  !
  REAL(DP), PARAMETER :: cnum(12) = [ 1._DP, -1._DP, 1._DP, -1._DP, 1._DP, &
    -691._DP, 1._DP, -3617._DP, 43867._DP, -174611._DP, 77683._DP, -236364091._DP ]
  REAL(DP), PARAMETER :: cdenom(12) = [ 12._DP, 120._DP, 252._DP, 240._DP, 132._DP, &
    32760._DP, 12._DP, 8160._DP, 14364._DP, 6600._DP, 276._DP, 65520._DP ]
  !* FIRST EXECUTABLE STATEMENT  DXPSI
  n = MAX(0,Ipsix-INT(A))
  b = n + A
  k1 = Ipsik - 1
  !
  !        SERIES EXPANSION FOR A > IPSIX USING IPSIK-1 TERMS.
  !
  c = 0._DP
  DO i = 1, k1
    k = Ipsik - i
    c = (c+cnum(k)/cdenom(k))/b**2
  END DO
  DXPSI = LOG(b) - (c+.5_DP/b)
  IF( n/=0 ) THEN
    b = 0._DP
    !
    !        RECURRENCE FOR A <= IPSIX.
    !
    DO m = 1, n
      b = b + 1._DP/(n-m+A)
    END DO
    DXPSI = DXPSI - b
  END IF
END FUNCTION DXPSI
