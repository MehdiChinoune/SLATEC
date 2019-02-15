!DECK MPMAXR
SUBROUTINE MPMAXR(X)
  IMPLICIT NONE
  INTEGER i, it, LUN, M, MXR
  !***BEGIN PROLOGUE  MPMAXR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPMAXR-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !  Sets X to the largest possible positive 'mp' number.
  !
  !  The argument X(*) is an INTEGER arrays of size 30.  See the comments
  !  in the routine MPBLAS for the reason for this choice.
  !
  !***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
  !***ROUTINES CALLED  MPCHK
  !***COMMON BLOCKS    MPCOM
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  !***END PROLOGUE  MPMAXR
  COMMON /MPCOM / B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*)
  !***FIRST EXECUTABLE STATEMENT  MPMAXR
  CALL MPCHK(1,4)
  it = B - 1
  ! SET FRACTION DIGITS TO B-1
  DO i = 1, T
    X(i+2) = it
  ENDDO
  ! SET SIGN AND EXPONENT
  X(1) = 1
  X(2) = M
END SUBROUTINE MPMAXR
