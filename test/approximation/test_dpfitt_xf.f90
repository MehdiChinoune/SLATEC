MODULE TEST_DPFITT_XF
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DPFITT()
    !> Quick check for DPOLFT, DPCOEF and DP1VLU.
    USE service, ONLY : eps_dp
    USE approximation, ONLY : DP1VLU, DPCOEF, DPOLFT
    !     .. Local Scalars ..
    REAL(DP) :: eps, rp, sveps, tol
    INTEGER :: i, m, maxord, ierp, ierr, nord, nordp
    !     .. Local Arrays ..
    REAL(DP) :: a(97), w(11), x(11), y(11), r(11)
    !* FIRST EXECUTABLE STATEMENT  DPFITT
    !
    !     Initialize variables for testing passage or failure of tests
    !
    tol = SQRT(eps_dp)
    m = 11
    DO i = 1, m
      x(i) = i - 6
      y(i) = x(i)**4
    END DO
    !
    !     Test DPOLFT
    !     Input EPS is negative - specified level
    !
    rp = 625._DP
    w(1) = -1._DP
    maxord = 5
    !
    !     Input EPS is positive
    !
    nordp = 4
    eps = 75._DP*eps_dp
    sveps = eps
    !
    !     Improper input
    !
    ierp = 2
    m = -2
    CALL DPOLFT(m,x,y,w,maxord,nord,eps,r,ierr,a)
    !
    RETURN
  END SUBROUTINE DPFITT
  !
END MODULE TEST_DPFITT_XF
!
PROGRAM MAIN
  USE TEST_DPFITT_XF
  IMPLICIT NONE
  !
  CALL DPFITT()
  !
END PROGRAM MAIN