MODULE TEST_DCOV_XF
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DNLS1Q()
    !> Quick check for DNLS1E, DNLS1 and DCOV.
    !***
    ! **Description:**
    !
    !   This subroutine performs a quick check on the subroutines DNLS1E
    !   (and DNLS1) and DCOV.

    USE service, ONLY : eps_dp
    USE approximation, ONLY : DCOV, DNLS1E
    !     .. Local Scalars ..
    REAL(DP) :: tol, tol2
    INTEGER :: iflag, info, iopt, ldfjac, m, n, nprint
    !     .. Local Arrays ..
    REAL(DP) :: fjac(10,2), fvec(10), wa(40), x(2)
    !* FIRST EXECUTABLE STATEMENT  DNLS1Q
    !
    m = 10
    n = 2
    ldfjac = 10
    nprint = -1
    iflag = 1
    tol = MAX(SQRT(40._DP*eps_dp),1.E-12_DP)
    tol2 = SQRT(tol)
    !
    !     Test improper input parameters.
    !
    iopt = 2
    x(1) = 3.0E-1_DP
    x(2) = 4.0E-1_DP
    m = 0
    CALL DCOV(DFCN2,iopt,m,n,x,fvec,fjac,ldfjac,info,wa(1),wa(n+1),wa(2*n+1),wa(3*n+1))
    !
    RETURN
  END SUBROUTINE DNLS1Q
  !** DFCN2
  PURE SUBROUTINE DFCN2(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
    !> Subsidiary to DNLS1Q.
    !***
    ! **Description:**
    !
    !   Subroutine to evaluate function and full Jacobian for test
    !   problem in quick check of DNLS1E.

    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: Iflag, Ldfjac, M, N
    !     .. Array Arguments ..
    REAL(DP), INTENT(IN) :: X(N)
    REAL(DP), INTENT(INOUT) :: Fvec(M)
    REAL(DP), INTENT(OUT) :: Fjac(:,:)
    !     .. Local Scalars ..
    REAL(DP) :: temp
    INTEGER :: i
    !* FIRST EXECUTABLE STATEMENT  DFCN2
    IF( Iflag==0 ) RETURN
    !
    !     Should we evaluate function or Jacobian?
    !
    IF( Iflag==1 ) THEN
      !
      !       Evaluate functions.
      !
      DO i = 1, M
        temp = i
        Fvec(i) = 2._DP + 2._DP*temp - EXP(temp*X(1)) - EXP(temp*X(2))
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
  END SUBROUTINE DFCN2
  !
END MODULE TEST_DCOV_XF
!
PROGRAM MAIN
  USE TEST_DCOV_XF
  IMPLICIT NONE
  !
  CALL DNLS1Q()
  !
END PROGRAM MAIN