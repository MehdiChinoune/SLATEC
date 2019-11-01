MODULE TEST_PFITQX_XF
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE PFITQX()
    !> Quick check for POLFIT, PCOEF and PVALUE.
    USE service, ONLY : eps_sp
    USE approximation, ONLY : PCOEF, POLFIT, PVALUE
    !     .. Local Scalars ..
    REAL(SP) :: eps, rp, sveps, tol
    INTEGER :: i, m, maxord, ierp, ierr, nord, nordp
    !     .. Local Arrays ..
    REAL(SP) :: a(97), w(11), x(11), y(11), r(11)
    !* FIRST EXECUTABLE STATEMENT  PFITQX
    !
    !     Initialize variables for testing passage or failure of tests
    !
    tol = SQRT(eps_sp)
    m = 11
    DO i = 1, m
      x(i) = i - 6
      y(i) = x(i)**4
    END DO
    !
    rp = 625._SP
    w(1) = -1._SP
    maxord = 5
    nordp = 4
    eps = 75._SP*eps_sp
    sveps = eps
    !
    !     Improper input
    !
    ierp = 2
    m = -2
    CALL POLFIT(m,x,y,w,maxord,nord,eps,r,ierr,a)
    !
    RETURN
  END SUBROUTINE PFITQX
  !
END MODULE TEST_PFITQX_XF
!
PROGRAM MAIN
  USE TEST_PFITQX_XF
  IMPLICIT NONE
  !
  CALL PFITQX()
  !
END PROGRAM MAIN