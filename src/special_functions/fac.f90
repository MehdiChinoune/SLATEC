!** FAC
REAL FUNCTION FAC(N)
  !>
  !***
  !  Compute the factorial function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C1
  !***
  ! **Type:**      SINGLE PRECISION (FAC-S, DFAC-D)
  !***
  ! **Keywords:**  FACTORIAL, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! FAC(N) evaluates the factorial function of N.  FAC is single
  ! precision.  N must be an integer between 0 and 25 inclusive.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  GAMLIM, R9LGMC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  REAL x, xmax, xmin
  INTEGER N
  REAL, PARAMETER :: facn(26) = [ 1.0E0, 1.0E0, 2.0E0, 6.0E0, 24.0E0, 120.0E0, &
    720.0E0, 5040.0E0, 40320.0E0, 362880.0E0, 3628800.0E0, 39916800.0E0, &
    479001600.0E0, 6227020800.0E0, 87178291200.0E0, 1307674368000.0E0, &
    20922789888000.0E0,     355687428096000.0E0,   6402373705728000.0E0, &
    .12164510040883200E18, .24329020081766400E19, .51090942171709440E20, &
    .11240007277776077E22, .25852016738884977E23, .62044840173323944E24, &
    .15511210043330986E26 ]
  REAL, PARAMETER :: sq2pil = 0.91893853320467274E0
  INTEGER :: nmax = 0
  !* FIRST EXECUTABLE STATEMENT  FAC
  IF ( nmax==0 ) THEN
    CALL GAMLIM(xmin,xmax)
    nmax = INT( xmax ) - 1
  END IF
  !
  IF ( N<0 ) CALL XERMSG('SLATEC','FAC',&
    'FACTORIAL OF NEGATIVE INTEGER UNDEFINED',1,2)
  !
  IF ( N<=25 ) FAC = facn(N+1)
  IF ( N<=25 ) RETURN
  !
  IF ( N>nmax ) CALL XERMSG('SLATEC','FAC','N SO BIG FACTORIAL(N) OVERFLOWS',2,2)
  !
  x = N + 1
  FAC = EXP((x-0.5)*LOG(x)-x+sq2pil+R9LGMC(x))
  !
END FUNCTION FAC
