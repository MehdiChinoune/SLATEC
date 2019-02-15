!DECK RGAUSS
FUNCTION RGAUSS(Xmean,Sd)
  IMPLICIT NONE
  INTEGER i
  REAL RAND, RGAUSS, Sd, Xmean
  !***BEGIN PROLOGUE  RGAUSS
  !***PURPOSE  Generate a normally distributed (Gaussian) random number.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  L6A14
  !***TYPE      SINGLE PRECISION (RGAUSS-S)
  !***KEYWORDS  FNLIB, GAUSSIAN, NORMAL, RANDOM NUMBER, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Generate a normally distributed random number, i.e., generate random
  ! numbers with a Gaussian distribution.  These random numbers are not
  ! exceptionally good -- especially in the tails of the distribution,
  ! but this implementation is simple and suitable for most applications.
  ! See R. W. Hamming, Numerical Methods for Scientists and Engineers,
  ! McGraw-Hill, 1962, pages 34 and 389.
  !
  !             Input Arguments --
  ! XMEAN  the mean of the Guassian distribution.
  ! SD     the standard deviation of the Guassian function
  !          EXP (-1/2 * (X-XMEAN)**2 / SD**2)
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  RAND
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   910819  Added EXTERNAL statement for RAND due to problem on IBM
  !           RS 6000.  (WRB)
  !***END PROLOGUE  RGAUSS
  EXTERNAL RAND
  !***FIRST EXECUTABLE STATEMENT  RGAUSS
  RGAUSS = -6.0
  DO i = 1, 12
    RGAUSS = RGAUSS + RAND(0.0)
  ENDDO
  !
  RGAUSS = Xmean + Sd*RGAUSS
  !
END FUNCTION RGAUSS
