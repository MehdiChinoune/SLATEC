!DECK ENORM
REAL FUNCTION ENORM(N,X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ENORM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SNLS1, SNLS1E, SNSQ and SNSQE
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (ENORM-S, DENORM-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Given an N-vector X, this function calculates the
  !     Euclidean norm of X.
  !
  !     The Euclidean norm is computed by accumulating the sum of
  !     squares in three different sums. The sums of squares for the
  !     small and large components are scaled so that no overflows
  !     occur. Non-destructive underflows are permitted. Underflows
  !     and overflows do not occur in the computation of the unscaled
  !     sum of squares for the intermediate components.
  !     The definitions of small, intermediate and large components
  !     depend on two constants, RDWARF and RGIANT. The main
  !     restrictions on these constants are that RDWARF**2 not
  !     underflow and RGIANT**2 not overflow. The constants
  !     given here are suitable for every known computer.
  !
  !     The function statement is
  !
  !       REAL FUNCTION ENORM(N,X)
  !
  !     where
  !
  !       N is a positive integer input variable.
  !
  !       X is an input array of length N.
  !
  !***SEE ALSO  SNLS1, SNLS1E, SNSQ, SNSQE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  ENORM
  INTEGER N
  REAL X(*)
  INTEGER i
  REAL agiant, floatn, one, rdwarf, rgiant, s1, s2, s3, xabs, &
    x1max, x3max, zero
  SAVE one, zero, rdwarf, rgiant
  DATA one, zero, rdwarf, rgiant/1.0E0, 0.0E0, 3.834E-20, 1.304E19/
  !***FIRST EXECUTABLE STATEMENT  ENORM
  s1 = zero
  s2 = zero
  s3 = zero
  x1max = zero
  x3max = zero
  floatn = N
  agiant = rgiant/floatn
  DO i = 1, N
    xabs = ABS(X(i))
    IF ( xabs>rdwarf.AND.xabs<agiant ) THEN
      !
      !           SUM FOR INTERMEDIATE COMPONENTS.
      !
      s2 = s2 + xabs**2
    ELSEIF ( xabs<=rdwarf ) THEN
      !
      !              SUM FOR SMALL COMPONENTS.
      !
      IF ( xabs<=x3max ) THEN
        IF ( xabs/=zero ) s3 = s3 + (xabs/x3max)**2
      ELSE
        s3 = one + s3*(x3max/xabs)**2
        x3max = xabs
      ENDIF
      !
      !              SUM FOR LARGE COMPONENTS.
      !
    ELSEIF ( xabs<=x1max ) THEN
      s1 = s1 + (xabs/x1max)**2
    ELSE
      s1 = one + s1*(x1max/xabs)**2
      x1max = xabs
    ENDIF
  ENDDO
  !
  !     CALCULATION OF NORM.
  !
  IF ( s1/=zero ) THEN
    ENORM = x1max*SQRT(s1+(s2/x1max)/x1max)
  ELSEIF ( s2==zero ) THEN
    ENORM = x3max*SQRT(s3)
  ELSE
    IF ( s2>=x3max ) ENORM = SQRT(s2*(one+(x3max/s2)*(x3max*s3)))
    IF ( s2<x3max ) ENORM = SQRT(x3max*((s2/x3max)+(x3max*s3)))
  ENDIF
  !
  !     LAST CARD OF FUNCTION ENORM.
  !
END FUNCTION ENORM
