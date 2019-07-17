!** MPMAXR
PURE SUBROUTINE MPMAXR(X)
  !> Subsidiary to DQDOTA and DQDOTI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      ALL (MPMAXR-A)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !  Sets X to the largest possible positive 'mp' number.
  !
  !  The argument X(*) is an INTEGER arrays of size 30.  See the comments
  !  in the routine MPBLAS for the reason for this choice.
  !
  !***
  ! **See also:**  DQDOTA, DQDOTI, MPBLAS
  !***
  ! **Routines called:**  MPCHK
  !***
  ! COMMON BLOCKS    MPCOM

  !* REVISION HISTORY  (YYMMDD)
  !   791001  DATE WRITTEN
  !   ??????  Modified for use with BLAS.  Blank COMMON changed to named
  !           COMMON.  R given dimension 12.
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
  USE MPCOM, ONLY : b_com, m_com, t_com, mxr_com

  INTEGER, INTENT(OUT) :: X(mxr_com)
  INTEGER :: i, it
  !* FIRST EXECUTABLE STATEMENT  MPMAXR
  it = b_com - 1
  ! SET FRACTION DIGITS TO B-1
  DO i = 1, t_com
    X(i+2) = it
  END DO
  ! SET SIGN AND EXPONENT
  X(1) = 1
  X(2) = m_com

END SUBROUTINE MPMAXR