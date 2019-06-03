!** DXPSI
REAL(DP) FUNCTION DXPSI(A,Ipsik,Ipsix)
  !>
  !  To compute values of the Psi function for DXLEGF.
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

  INTEGER i, Ipsik, Ipsix, k, k1, m, n
  REAL(DP) :: A, b, c
  !
  !        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
  !        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
  !        NUMBER.
  !
  REAL(DP), PARAMETER :: cnum(12) = [ 1.D0, -1.D0, 1.D0, -1.D0, 1.D0, -691.D0, 1.D0, -3617.D0, &
    43867.D0, -174611.D0, 77683.D0, -236364091.D0 ]
  REAL(DP), PARAMETER :: cdenom(12) = [ 12.D0, 120.D0, 252.D0, 240.D0, 132.D0, 32760.D0, &
    12.D0, 8160.D0, 14364.D0, 6600.D0, 276.D0, 65520.D0 ]
  !* FIRST EXECUTABLE STATEMENT  DXPSI
  n = MAX(0,Ipsix-INT(A))
  b = n + A
  k1 = Ipsik - 1
  !
  !        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS.
  !
  c = 0.D0
  DO i = 1, k1
    k = Ipsik - i
    c = (c+cnum(k)/cdenom(k))/b**2
  END DO
  DXPSI = LOG(b) - (c+.5D0/b)
  IF ( n/=0 ) THEN
    b = 0.D0
    !
    !        RECURRENCE FOR A .LE. IPSIX.
    !
    DO m = 1, n
      b = b + 1.D0/(n-m+A)
    END DO
    DXPSI = DXPSI - b
  END IF
END FUNCTION DXPSI
