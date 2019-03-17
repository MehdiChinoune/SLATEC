!DECK DCSEVL
REAL(8) FUNCTION DCSEVL(X,Cs,N)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DCSEVL
  !***PURPOSE  Evaluate a Chebyshev series.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C3A2
  !***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
  !***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !  Evaluate the N-term Chebyshev series CS at X.  Adapted from
  !  a method presented in the paper by Broucke referenced below.
  !
  !       Input Arguments --
  !  X    value at which the series is to be evaluated.
  !  CS   array of N terms of a Chebyshev series.  In evaluating
  !       CS, only half the first coefficient is summed.
  !  N    number of terms in array CS.
  !
  !***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
  !                 Chebyshev series, Algorithm 446, Communications of
  !                 the A.C.M. 16, (1973) pp. 254-256.
  !               L. Fox and I. B. Parker, Chebyshev Polynomials in
  !                 Numerical Analysis, Oxford University Press, 1968,
  !                 page 56.
  !***ROUTINES CALLED  D1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900329  Prologued revised extensively and code rewritten to allow
  !           X to be slightly outside interval (-1,+1).  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DCSEVL
  INTEGER i, N, ni
  REAL(8) :: b0, b1, b2, Cs(*), onepl, twox, X, D1MACH
  LOGICAL first
  SAVE first, onepl
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DCSEVL
  IF ( first ) onepl = 1.0D0 + D1MACH(4)
  first = .FALSE.
  IF ( N<1 ) CALL XERMSG('SLATEC','DCSEVL','NUMBER OF TERMS .LE. 0',2,2)
  IF ( N>1000 ) CALL XERMSG('SLATEC','DCSEVL','NUMBER OF TERMS .GT. 1000',3,&
    2)
  IF ( ABS(X)>onepl ) CALL XERMSG('SLATEC','DCSEVL',&
    'X OUTSIDE THE INTERVAL (-1,+1)',1,1)
  !
  b1 = 0.0D0
  b0 = 0.0D0
  twox = 2.0D0*X
  DO i = 1, N
    b2 = b1
    b1 = b0
    ni = N + 1 - i
    b0 = twox*b1 - b2 + Cs(ni)
  ENDDO
  !
  DCSEVL = 0.5D0*(b0-b2)
  !
END FUNCTION DCSEVL