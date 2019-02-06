!*==BINOM.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK BINOM
      FUNCTION BINOM(N,M)
      IMPLICIT NONE
!*--BINOM5
!*** Start of declarations inserted by SPAG
      REAL ALNREL , bilnmx , BINOM , corr , fintmx , R1MACH , R9LGMC , sq2pil , 
     &     xk , xn , xnk
      INTEGER i , k , M , N
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  BINOM
!***PURPOSE  Compute the binomial coefficients.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1
!***TYPE      SINGLE PRECISION (BINOM-S, DBINOM-D)
!***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! BINOM(N,M) calculates the binomial coefficient (N!)/((M!)*(N-M)!).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALNREL, R1MACH, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  BINOM
      LOGICAL first
      SAVE sq2pil , bilnmx , fintmx , first
      DATA sq2pil/0.91893853320467274E0/
      DATA first/.TRUE./
!***FIRST EXECUTABLE STATEMENT  BINOM
      IF ( first ) THEN
        bilnmx = LOG(R1MACH(2))
        fintmx = 0.9/R1MACH(3)
      ENDIF
      first = .FALSE.
!
      IF ( N<0.OR.M<0 ) CALL XERMSG('SLATEC','BINOM','N OR M LT ZERO',1,2)
      IF ( N<M ) CALL XERMSG('SLATEC','BINOM','N LT M',2,2)
!
      k = MIN(M,N-M)
      IF ( k<=20 ) THEN
        IF ( k*LOG(AMAX0(N,1))<=bilnmx ) THEN
!
          BINOM = 1.
          IF ( k==0 ) RETURN
!
          DO i = 1 , k
            BINOM = BINOM*REAL(N-i+1)/i
          ENDDO
!
          IF ( BINOM<fintmx ) BINOM = AINT(BINOM+0.5)
          RETURN
        ENDIF
      ENDIF
!
! IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
      IF ( k<9 ) CALL XERMSG('SLATEC','BINOM',
     &                       'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG',3,2)
!
      xn = N + 1
      xk = k + 1
      xnk = N - k + 1
!
      corr = R9LGMC(xn) - R9LGMC(xk) - R9LGMC(xnk)
      BINOM = xk*LOG(xnk/xk) - xn*ALNREL(-(xk-1.)/xn) - 0.5*LOG(xn*xnk/xk)
     &        + 1.0 - sq2pil + corr
!
      IF ( BINOM>bilnmx ) CALL XERMSG('SLATEC','BINOM',
     &                             'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG'
     &                             ,3,2)
!
      BINOM = EXP(BINOM)
      IF ( BINOM<fintmx ) BINOM = AINT(BINOM+0.5)
!
      END FUNCTION BINOM
