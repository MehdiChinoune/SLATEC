MODULE TEST_DPSORT_XF_1
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DSRTQC()
    !> Quick check for SLATEC routines DSORT, DPSORT, DPPERM
    USE data_handling, ONLY : DPSORT
    !
    INTEGER, PARAMETER :: N = 9
    !
    REAL(DP) :: y(N)
    INTEGER :: iy(N), ier, nn, kkflag
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
    CALL DPSORT(y,N,iy,1,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = -1
    kkflag = 1
    CALL DPSORT(y,nn,iy,kkflag,ier)
    !
    RETURN
  END SUBROUTINE DSRTQC
  !
END MODULE TEST_DPSORT_XF_1
!
PROGRAM MAIN
  USE TEST_DPSORT_XF_1
  IMPLICIT NONE
  !
  CALL DSRTQC()
  !
END PROGRAM MAIN