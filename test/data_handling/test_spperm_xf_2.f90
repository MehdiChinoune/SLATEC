MODULE TEST_SPPERM_XF_2
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE SSRTQC()
    !> Quick check for SLATEC routines SSORT, SPSORT, SPPERM
    USE data_handling, ONLY : SPPERM
    !
    INTEGER, PARAMETER :: N = 9
    !
    REAL(SP) :: y(N)
    INTEGER :: iy(N), i, ier, nn
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !
    REAL(SP), PARAMETER :: x(N) = [ 36._SP, 54._SP, -1._SP, 29._SP, 1._SP, &
      80._SP, 98._SP, 99._SP, 55._SP ]
    !* FIRST EXECUTABLE STATEMENT  SSRTQC
    !
    !        ... SETUP PROBLEM
    !
    y = x
    iy = [ ( i, i = 1, N ) ]
    CALL SPPERM(y,N,iy,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = 1
    iy(1) = 5
    CALL SPPERM(y,nn,iy,ier)
    !
    RETURN
  END SUBROUTINE SSRTQC
  !
END MODULE TEST_SPPERM_XF_2
!
PROGRAM MAIN
  USE TEST_SPPERM_XF_2
  IMPLICIT NONE
  !
  CALL SSRTQC()
  !
END PROGRAM MAIN