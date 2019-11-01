MODULE TEST_HPPERM_XF_2
  IMPLICIT NONE

CONTAINS
  SUBROUTINE HSRTQC()
    !> Quick check for SLATEC routine HPSORT, HPPERM
    USE data_handling, ONLY : HPPERM
    !
    INTEGER, PARAMETER :: N = 9
    !
    CHARACTER(2) :: y(N), work
    INTEGER :: iy(N), i, ier, nn
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !
    CHARACTER(2), PARAMETER :: x(N) = [ 'AC', 'AZ', 'AD', 'AA', 'AB', 'ZZ', &
      'ZA', 'ZX', 'ZY' ]
    !* FIRST EXECUTABLE STATEMENT  HSRTQC
    !
    !        ... SETUP PROBLEM
    !
    y = x
    iy = [ ( i, i = 1, N ) ]
    CALL HPPERM(y,N,iy,work,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = 1
    iy(1) = 5
    CALL HPPERM(y,nn,iy,work,ier)
    !
    RETURN
  END SUBROUTINE HSRTQC
  !
END MODULE TEST_HPPERM_XF_2
!
PROGRAM MAIN
  USE TEST_HPPERM_XF_2
  IMPLICIT NONE
  !
  CALL HSRTQC()
  !
END PROGRAM MAIN