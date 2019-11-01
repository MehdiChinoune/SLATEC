MODULE TEST_SCOV_XF
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE SNLS1Q()
    !> Quick check for SNLS1E, SNLS1 and SCOV.
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutines SNLS1E
    !   (and SNLS1) and SCOV.

    USE service, ONLY : eps_sp
    USE approximation, ONLY : SCOV
    !     .. Local Scalars ..
    REAL(SP) :: tol
    INTEGER :: info, iopt, ldfjac, lwa, m, n, nprint
    !     .. Local Arrays ..
    REAL(SP) :: fjac(10,2), fvec(10), wa(40), x(2)
    !* FIRST EXECUTABLE STATEMENT  SNLS1Q
    m = 10
    n = 2
    lwa = 40
    ldfjac = 10
    nprint = -1
    tol = SQRT(40._SP*eps_sp)
    !
    !     Test improper input parameters.
    !
    lwa = 35
    iopt = 2
    x(1) = 3.0E-1_SP
    x(2) = 4.0E-1_SP
    m = 0
    CALL SCOV(FCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),wa(3*n+1))
    !
    RETURN
  END SUBROUTINE SNLS1Q
  !** FCN2
  PURE SUBROUTINE FCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    !> Subsidiary to SNLS1Q.
    !***
    ! **Description:**
    !
    !   Subroutine to evaluate function and full Jacobian for test
    !   problem in quick check of SNLS1E.

    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL(SP), INTENT(IN) :: X(N)
    REAL(SP), INTENT(INOUT) :: Fvec(M)
    REAL(SP), INTENT(OUT) :: Fjac(:,:)
    !     .. Local Scalars ..
    REAL(SP) :: temp
    INTEGER :: i
    !* FIRST EXECUTABLE STATEMENT  FCN2
    IF( Iflag==0 ) RETURN
    !
    !     Should we evaluate functions or Jacobian?
    !
    IF( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = 2._SP + 2._SP*temp - EXP(temp*X(1)) - EXP(temp*X(2))
      END DO
    ELSE
      !
      !       Evaluate Jacobian.
      !
      IF( Iflag/=2 ) RETURN
      DO i = 1, M
        temp = i
        Fjac(i,1) = -temp*EXP(temp*X(1))
        Fjac(i,2) = -temp*EXP(temp*X(2))
      END DO
    END IF
  END SUBROUTINE FCN2
  !
END MODULE TEST_SCOV_XF
!
PROGRAM MAIN
  USE TEST_SCOV_XF
  IMPLICIT NONE
  !
  CALL SNLS1Q()
  !
END PROGRAM MAIN