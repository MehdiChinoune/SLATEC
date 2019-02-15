!DECK CBEG
COMPLEX FUNCTION CBEG(Reset)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CBEG
  !***SUBSIDIARY
  !***PURPOSE  Generate random numbers.
  !***LIBRARY   SLATEC (BLAS)
  !***AUTHOR  Du Croz, J. (NAG)
  !           Hanson, R. J. (SNLA)
  !***DESCRIPTION
  !
  !  Generates random numbers uniformly distributed between -0.5 and 0.5.
  !
  !  Auxiliary routine for test program for Level 2 Blas.
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   870810  DATE WRITTEN
  !   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
  !***END PROLOGUE  CBEG
  !     .. Scalar Arguments ..
  LOGICAL Reset
  !     .. Local Scalars ..
  INTEGER i, ic, j, mi, mj
  !     .. Save statement ..
  SAVE i, ic, j, mi, mj
  !     .. Intrinsic Functions ..
  INTRINSIC CMPLX
  !***FIRST EXECUTABLE STATEMENT  CBEG
  IF ( Reset ) THEN
    !        Initialize local variables.
    mi = 891
    mj = 457
    i = 7
    j = 7
    ic = 0
    Reset = .FALSE.
  ENDIF
  !
  !     The sequence of values of I or J is bounded between 1 and 999.
  !     If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
  !     If initial I or J = 4 or 8, the period will be 25.
  !     If initial I or J = 5, the period will be 10.
  !     IC is used to break up the period by skipping 1 value of I or J
  !     in 6.
  !
  ic = ic + 1
  DO
    i = i*mi
    j = j*mj
    i = i - 1000*(i/1000)
    j = j - 1000*(j/1000)
    IF ( ic>=5 ) THEN
      ic = 0
      CYCLE
    ENDIF
    CBEG = CMPLX((i-500)/1001.0,(j-500)/1001.0)
    EXIT
  ENDDO
  !
  !     End of CBEG.
  !
END FUNCTION CBEG
