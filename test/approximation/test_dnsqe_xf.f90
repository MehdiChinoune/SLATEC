MODULE TEST_DNSQE_XF
  USE service, ONLY : DP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE DNSQQK()
    !> Quick check for DNSQE and DNSQ.
    !***
    ! **Description:**
    !
    !  This subroutine performs a quick check on the subroutine DNSQE (and DNSQ).
    USE service, ONLY : eps_dp
    USE approximation, ONLY : DNSQE
    !     .. Local Scalars ..
    REAL(DP) :: tol
    INTEGER :: info, iopt, lwa, n, nprint
    !     .. Local Arrays ..
    REAL(DP) :: fvec(2), wa(19), x(2)
    !* FIRST EXECUTABLE STATEMENT  DNSQQK
    n = 2
    nprint = -1
    tol = SQRT(eps_dp)
    !
    !     Test improper input parameters.
    !
    lwa = 15
    iopt = 1
    x(1) = -1.2_DP
    x(2) = 1._DP
    CALL DNSQE(DQFCN2,DQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
    !
    RETURN
  END SUBROUTINE DNSQQK
  !** DQFCN2
  PURE SUBROUTINE DQFCN2(N,X,Fvec,Iflag)
    !> Evaluate function used in DNSQE.
    !***
    ! **Description:**
    !
    !   Subroutine which evaluates the function for test program
    !   used in quick check of DNSQE.

    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: Iflag, N
    !     .. Array Arguments ..
    REAL(DP), INTENT(IN) :: X(N)
    REAL(DP), INTENT(OUT) :: Fvec(N)
    !* FIRST EXECUTABLE STATEMENT  DQFCN2
    Fvec(1) = 1._DP - X(1)
    Fvec(2) = 10._DP*(X(2)-X(1)**2)
  END SUBROUTINE DQFCN2
  !** DQJAC2
  PURE SUBROUTINE DQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
    !>
    !***
    ! **Description:**
    !
    !  Subroutine to evaluate the full jacobian for test problem used in quick
    !  check of DNSQE.

    INTEGER, INTENT(IN) :: Iflag, Ldfjac, N
    REAL(DP), INTENT(IN) :: Fvec(N), X(N)
    REAL(DP), INTENT(OUT) :: Fjac(Ldfjac,N)
    !* FIRST EXECUTABLE STATEMENT  DQJAC2
    Fjac(1,1) = -1._DP
    Fjac(1,2) = 0._DP
    Fjac(2,1) = -20._DP*X(1)
    Fjac(2,2) = 10._DP
  END SUBROUTINE DQJAC2
  !
END MODULE TEST_DNSQE_XF
!
PROGRAM MAIN
  USE TEST_DNSQE_XF
  IMPLICIT NONE
  !
  CALL DNSQQK()
  !
END PROGRAM MAIN