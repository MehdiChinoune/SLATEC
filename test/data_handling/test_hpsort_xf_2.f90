MODULE TEST_HPSORT_XF_2
  IMPLICIT NONE

CONTAINS
  SUBROUTINE HSRTQC()
    !> Quick check for SLATEC routine HPSORT, HPPERM
    USE data_handling, ONLY : HPSORT
    !
    INTEGER, PARAMETER :: N = 9
    !
    CHARACTER(2) :: y(N), work
    INTEGER :: iy(N), ier, strbeg, strend, nn, kkflag
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
    CALL HPSORT(y,N,1,2,iy,1,work,ier)
    !
    !     ... TEST ERROR MESSAGES
    !
    nn = 1
    strbeg = 1
    strend = 2
    kkflag = 0
    CALL HPSORT(y,nn,strbeg,strend,iy,kkflag,work,ier)
    !
    RETURN
  END SUBROUTINE HSRTQC
  !
END MODULE TEST_HPSORT_XF_2
!
PROGRAM MAIN
  USE TEST_HPSORT_XF_2
  IMPLICIT NONE
  !
  CALL HSRTQC()
  !
END PROGRAM MAIN