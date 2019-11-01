MODULE TEST_IPSORT_XF_1
  IMPLICIT NONE

CONTAINS
  SUBROUTINE ISRTQC()
    !> Quick check for SLATEC routines ISORT, IPSORT, IPPERM
    USE data_handling, ONLY : IPSORT
    !
    INTEGER, PARAMETER :: N = 9
    !
    INTEGER :: y(N), iy(N)
    INTEGER :: ier, nn, kkflag
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
    !                            CHECK IPSORT
    !     -------------------------------------------------------------
    y = x
    CALL IPSORT(y,N,iy,1,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = -1
    kkflag = 1
    CALL IPSORT(y,nn,iy,kkflag,ier)
    !
  END SUBROUTINE ISRTQC
  !
END MODULE TEST_IPSORT_XF_1
!
PROGRAM MAIN
  USE TEST_IPSORT_XF_1
  IMPLICIT NONE
  !
  CALL ISRTQC()
  !
END PROGRAM MAIN