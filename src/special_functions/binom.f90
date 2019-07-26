!** BINOM
REAL(SP) ELEMENTAL FUNCTION BINOM(N,M)
  !> Compute the binomial coefficients.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C1
  !***
  ! **Type:**      SINGLE PRECISION (BINOM-S, DBINOM-D)
  !***
  ! **Keywords:**  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BINOM(N,M) calculates the binomial coefficient (N!)/((M!)*(N-M)!).
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  ALNREL, R1MACH, R9LGMC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  USE service, ONLY : huge_sp, eps_2_sp
  !
  INTEGER, INTENT(IN) :: M, N
  !
  INTEGER :: i, k
  REAL(SP) :: corr, xk, xn, xnk
  REAL(SP), PARAMETER :: bilnmx = LOG(huge_sp), fintmx = 0.9_SP/eps_2_sp
  REAL(SP), PARAMETER :: sq2pil = 0.91893853320467274_SP
  !* FIRST EXECUTABLE STATEMENT  BINOM
  !
  IF( N<0 .OR. M<0 ) THEN
    ERROR STOP 'BINOM : N OR M < 0'
  ELSEIF( N<M ) THEN
    ERROR STOP 'BINOM : N < M'
  END IF
  !
  k = MIN(M,N-M)
  IF( k<=20 ) THEN
    IF( k*LOG(AMAX0(N,1))<=bilnmx ) THEN
      !
      BINOM = 1._SP
      IF( k==0 ) RETURN
      !
      DO i = 1, k
        BINOM = BINOM*REAL(N-i+1,SP)/i
      END DO
      !
      IF( BINOM<fintmx ) BINOM = AINT(BINOM+0.5_SP)
      RETURN
    END IF
  END IF
  !
  ! IF K<9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
  IF( k<9 ) ERROR STOP 'BINOM : RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG'
  !
  xn = N + 1
  xk = k + 1
  xnk = N - k + 1
  !
  corr = R9LGMC(xn) - R9LGMC(xk) - R9LGMC(xnk)
  BINOM = xk*LOG(xnk/xk) - xn*ALNREL(-(xk-1._SP)/xn) - 0.5_SP*LOG(xn*xnk/xk)&
    + 1._SP - sq2pil + corr
  !
  IF( BINOM>bilnmx ) ERROR STOP 'BINOM : RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG'
  !
  BINOM = EXP(BINOM)
  IF( BINOM<fintmx ) BINOM = AINT(BINOM+0.5_SP)
  !
END FUNCTION BINOM