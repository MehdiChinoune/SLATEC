!*==FAC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK FAC
FUNCTION FAC(N)
  IMPLICIT NONE
  !*--FAC5
  !*** Start of declarations inserted by SPAG
  REAL FAC , facn , R9LGMC , sq2pil , x , xmax , xmin
  INTEGER N , nmax
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  FAC
  !***PURPOSE  Compute the factorial function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C1
  !***TYPE      SINGLE PRECISION (FAC-S, DFAC-D)
  !***KEYWORDS  FACTORIAL, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! FAC(N) evaluates the factorial function of N.  FAC is single
  ! precision.  N must be an integer between 0 and 25 inclusive.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  GAMLIM, R9LGMC, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  FAC
  DIMENSION facn(26)
  SAVE facn , sq2pil , nmax
  DATA facn(1)/1.0E0/
  DATA facn(2)/1.0E0/
  DATA facn(3)/2.0E0/
  DATA facn(4)/6.0E0/
  DATA facn(5)/24.0E0/
  DATA facn(6)/120.0E0/
  DATA facn(7)/720.0E0/
  DATA facn(8)/5040.0E0/
  DATA facn(9)/40320.0E0/
  DATA facn(10)/362880.0E0/
  DATA facn(11)/3628800.0E0/
  DATA facn(12)/39916800.0E0/
  DATA facn(13)/479001600.0E0/
  DATA facn(14)/6227020800.0E0/
  DATA facn(15)/87178291200.0E0/
  DATA facn(16)/1307674368000.0E0/
  DATA facn(17)/20922789888000.0E0/
  DATA facn(18)/355687428096000.0E0/
  DATA facn(19)/6402373705728000.0E0/
  DATA facn(20)/.12164510040883200E18/
  DATA facn(21)/.24329020081766400E19/
  DATA facn(22)/.51090942171709440E20/
  DATA facn(23)/.11240007277776077E22/
  DATA facn(24)/.25852016738884977E23/
  DATA facn(25)/.62044840173323944E24/
  DATA facn(26)/.15511210043330986E26/
  DATA sq2pil/0.91893853320467274E0/
  DATA nmax/0/
  !***FIRST EXECUTABLE STATEMENT  FAC
  IF ( nmax==0 ) THEN
    CALL GAMLIM(xmin,xmax)
    nmax = xmax - 1.
  ENDIF
  !
  IF ( N<0 ) CALL XERMSG('SLATEC','FAC',&
    'FACTORIAL OF NEGATIVE INTEGER UNDEFINED',1,2)
  !
  IF ( N<=25 ) FAC = facn(N+1)
  IF ( N<=25 ) RETURN
  !
  IF ( N>nmax ) CALL XERMSG('SLATEC','FAC','N SO BIG FACTORIAL(N) OVERFLOWS'&
    ,2,2)
  !
  x = N + 1
  FAC = EXP((x-0.5)*LOG(x)-x+sq2pil+R9LGMC(x))
  !
END FUNCTION FAC
