!** XPSI
REAL FUNCTION XPSI(A,Ipsik,Ipsix)
  IMPLICIT NONE
  !>
  !***
  !  To compute values of the Psi function for XLEGF.
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
  !           Corrected order of sections in prologue and added TYPE
  !           section.  (WRB)
  !   920127  Revised PURPOSE section of prologue.  (DWL)

  INTEGER i, Ipsik, Ipsix, k, k1, m, n
  REAL A, b, c
  !
  !        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
  !        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
  !        NUMBER.
  !
  REAL, PARAMETER :: cnum(12) = [ 1., -1., 1., -1., 1., -691., 1., -3617., 43867., &
    -174611., 77683., -236364091. ]
  REAL, PARAMETER :: cdenom(12) = [ 12., 120., 252., 240., 132., 32760., 12., 8160., &
    14364., 6600., 276., 65520. ]
  !* FIRST EXECUTABLE STATEMENT  XPSI
  n = MAX(0,Ipsix-INT(A))
  b = n + A
  k1 = Ipsik - 1
  !
  !        SERIES EXPANSION FOR A .GT. IPSIX USING IPSIK-1 TERMS.
  !
  c = 0.
  DO i = 1, k1
    k = Ipsik - i
    c = (c+cnum(k)/cdenom(k))/b**2
  ENDDO
  XPSI = LOG(b) - (c+.5/b)
  IF ( n/=0 ) THEN
    b = 0.
    !
    !        RECURRENCE FOR A .LE. IPSIX.
    !
    DO m = 1, n
      b = b + 1./(n-m+A)
    ENDDO
    XPSI = XPSI - b
  ENDIF
END FUNCTION XPSI
