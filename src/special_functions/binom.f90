!** BINOM
REAL FUNCTION BINOM(N,M)
  !>
  !  Compute the binomial coefficients.
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
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL corr, xk, xn, xnk
  INTEGER i, k, M, N
  REAL, SAVE :: bilnmx, fintmx
  REAL, PARAMETER :: sq2pil = 0.91893853320467274E0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BINOM
  IF ( first ) THEN
    bilnmx = LOG(R1MACH(2))
    fintmx = 0.9/R1MACH(3)
    first = .FALSE.
  END IF
  !
  IF ( N<0.OR.M<0 ) CALL XERMSG('BINOM','N OR M LT ZERO',1,2)
  IF ( N<M ) CALL XERMSG('BINOM','N LT M',2,2)
  !
  k = MIN(M,N-M)
  IF ( k<=20 ) THEN
    IF ( k*LOG(AMAX0(N,1))<=bilnmx ) THEN
      !
      BINOM = 1.
      IF ( k==0 ) RETURN
      !
      DO i = 1, k
        BINOM = BINOM*REAL(N-i+1)/i
      END DO
      !
      IF ( BINOM<fintmx ) BINOM = AINT(BINOM+0.5)
      RETURN
    END IF
  END IF
  !
  ! IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
  IF ( k<9 ) CALL XERMSG('BINOM',&
    'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG',3,2)
  !
  xn = N + 1
  xk = k + 1
  xnk = N - k + 1
  !
  corr = R9LGMC(xn) - R9LGMC(xk) - R9LGMC(xnk)
  BINOM = xk*LOG(xnk/xk) - xn*ALNREL(-(xk-1.)/xn) - 0.5*LOG(xn*xnk/xk)&
    + 1.0 - sq2pil + corr
  !
  IF ( BINOM>bilnmx ) CALL XERMSG('BINOM',&
    'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG',3,2)
  !
  BINOM = EXP(BINOM)
  IF ( BINOM<fintmx ) BINOM = AINT(BINOM+0.5)
  !
END FUNCTION BINOM
