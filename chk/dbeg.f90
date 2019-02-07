!*==DBEG.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DBEG
REAL(8) FUNCTION DBEG(Reset)
  IMPLICIT NONE
  !*--DBEG5
  !***BEGIN PROLOGUE  DBEG
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
  !***END PROLOGUE  DBEG
  !     .. Scalar Arguments ..
  LOGICAL Reset
  !     .. Local Scalars ..
  INTEGER i, ic, mi
  !     .. Save statement ..
  SAVE i, ic, mi
  !     .. Intrinsic Functions ..
  INTRINSIC REAL
  !***FIRST EXECUTABLE STATEMENT  DBEG
  IF ( Reset ) THEN
    !        Initialize local variables.
    mi = 891
    i = 7
    ic = 0
    Reset = .FALSE.
  ENDIF
  !
  !     The sequence of values of I is bounded between 1 and 999.
  !     If initial I = 1,2,3,6,7 or 9, the period will be 50.
  !     If initial I = 4 or 8, the period will be 25.
  !     If initial I = 5, the period will be 10.
  !     IC is used to break up the period by skipping 1 value of I in 6.
  !
  ic = ic + 1
  DO
    i = i*mi
    i = i - 1000*(i/1000)
    IF ( ic>=5 ) THEN
      ic = 0
      CYCLE
    ENDIF
    DBEG = REAL(i-500, 8)/1001.0D0
    EXIT
  ENDDO
  !
  !     End of DBEG.
  !
END FUNCTION DBEG
