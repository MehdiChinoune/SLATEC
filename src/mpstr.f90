!*==MPSTR.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK MPSTR
SUBROUTINE MPSTR(X,Y)
  IMPLICIT NONE
  !*--MPSTR5
  !*** Start of declarations inserted by SPAG
  INTEGER i, LUN, M, MXR
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  MPSTR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DQDOTA and DQDOTI
  !***LIBRARY   SLATEC
  !***TYPE      ALL (MPSTR-A)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !  Sets Y = X for 'mp' X and Y.
  !
  !  The arguments X(*) and Y(*) are INTEGER arrays of size 30.  See the
  !  comments in the routine MPBLAS for the reason for this choice.
  !
  !***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
  !***ROUTINES CALLED  (NONE)
  !***COMMON BLOCKS    MPCOM
  !***REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  !***END PROLOGUE  MPSTR
  COMMON /MPCOM / B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Y(*)
  !***FIRST EXECUTABLE STATEMENT  MPSTR
  DO i = 1, T + 2
    Y(i) = X(i)
  ENDDO
END SUBROUTINE MPSTR
