MODULE TEST_IPPERM_XF_2
  IMPLICIT NONE

CONTAINS
  SUBROUTINE ISRTQC()
    !> Quick check for SLATEC routines ISORT, IPSORT, IPPERM
    USE data_handling, ONLY : IPPERM
    !
    INTEGER, PARAMETER :: N = 9
    !
    INTEGER :: y(N), iy(N)
    INTEGER :: i, ier, nn
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !
    INTEGER, PARAMETER :: x(N) = [ 36, 54, -1, 29, 1, 80, 98, 99, 55 ]
    !* FIRST EXECUTABLE STATEMENT  ISRTQC
    !     -------------------------------------------------------------
    !                            CHECK IPPERM
    !     -------------------------------------------------------------
    y = x
    iy = [ ( i, i = 1, N ) ]
    CALL IPPERM(y,N,iy,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = 1
    iy(1) = 5
    CALL IPPERM(y,nn,iy,ier)
    !
  END SUBROUTINE ISRTQC
  !
END MODULE TEST_IPPERM_XF_2
!
PROGRAM MAIN
  USE TEST_IPPERM_XF_2
  IMPLICIT NONE
  !
  CALL ISRTQC()
  !
END PROGRAM MAIN