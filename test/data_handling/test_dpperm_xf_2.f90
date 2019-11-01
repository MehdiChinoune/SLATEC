MODULE TEST_DPPERM_XF_2
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DSRTQC()
    !> Quick check for SLATEC routines DSORT, DPSORT, DPPERM
    USE data_handling, ONLY : DPPERM
    !
    INTEGER, PARAMETER :: N = 9
    !
    REAL(DP) :: y(N)
    INTEGER :: iy(N), i, ier, nn
    !
    !     ---------
    !     TEST DATA
    !     ---------
    !
    !         X   = TEST VECTOR
    !
    REAL(DP), PARAMETER ::  x(N) = [ 36._DP, 54._DP, -1._DP, 29._DP, 1._DP, &
      80._DP, 98._DP, 99._DP, 55._DP ]
    !* FIRST EXECUTABLE STATEMENT  DSRTQC
    !
    !        ... SETUP PROBLEM
    !
    y = x
    iy = [ ( i, i = 1, N ) ]
    CALL DPPERM(y,N,iy,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = 1
    iy(1) = 5
    CALL DPPERM(y,nn,iy,ier)
    !
    RETURN
  END SUBROUTINE DSRTQC
  !
END MODULE TEST_DPPERM_XF_2
!
PROGRAM MAIN
  USE TEST_DPPERM_XF_2
  IMPLICIT NONE
  !
  CALL DSRTQC()
  !
END PROGRAM MAIN