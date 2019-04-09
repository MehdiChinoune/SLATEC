!** MPSTR
SUBROUTINE MPSTR(X,Y)
  USE MPCOM
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPSTR-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Sets Y = X for 'mp' X and Y.
  !
  !  The arguments X(*) and Y(*) are INTEGER arrays of size 30.  See the
  !  comments in the routine MPBLAS for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  (NONE)
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)

  INTEGER i, X(*), Y(*)
  !* FIRST EXECUTABLE STATEMENT  MPSTR
  DO i = 1, T + 2
    Y(i) = X(i)
  END DO
END SUBROUTINE MPSTR
