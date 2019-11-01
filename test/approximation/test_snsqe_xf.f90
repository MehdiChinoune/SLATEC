MODULE TEST_SNSQE_XF
  USE service, ONLY : SP
  IMPLICIT NONE

CONTAINS
  SUBROUTINE SNSQQK()
    !> Quick check for SNSQE and SNSQ.
    !***
    ! **Description:**
    !
    !  This subroutine performs a quick check on the subroutine SNSQE (and SNSQ).
    USE service, ONLY : eps_sp
    USE approximation, ONLY : SNSQE
    !     .. Local Scalars ..
    REAL(SP) :: tol
    INTEGER :: info, iopt, lwa, n, nprint
    !     .. Local Arrays ..
    REAL(SP) :: fvec(2), wa(19), x(2)
    !* FIRST EXECUTABLE STATEMENT  SNSQQK
    n = 2
    nprint = -1
    tol = SQRT(eps_sp)
    !
    !     Test improper input parameters.
    !
    lwa = 15
    iopt = 1
    x(1) = -1.2_SP
    x(2) = 1._SP
    CALL SNSQE(SQFCN2,SQJAC2,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
    !
    RETURN
  END SUBROUTINE SNSQQK
  !** SQFCN2
  PURE SUBROUTINE SQFCN2(N,X,Fvec,Iflag)
    !> Evaluate function used in SNSQE.
    !***
    ! **Description:**
    !
    !  Subroutine which evaluates the function for test program used in quick
    !  check of SNSQE.

    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: Iflag, N
    !     .. Array Arguments ..
    REAL(SP), INTENT(IN) :: X(N)
    REAL(SP), INTENT(OUT) :: Fvec(N)
    !* FIRST EXECUTABLE STATEMENT  SQFCN2
    Fvec(1) = 1._SP - X(1)
    Fvec(2) = 10._SP*(X(2)-X(1)**2)
  END SUBROUTINE SQFCN2
  !** SQJAC2
  PURE SUBROUTINE SQJAC2(N,X,Fvec,Fjac,Ldfjac,Iflag)
    !> Evaluate full Jacobian for SNSQE test.
    !***
    ! **Description:**
    !
    !  Subroutine to evaluate the full jacobian for test problem used in quick
    !  check of SNSQE.

    INTEGER, INTENT(IN) :: Iflag, Ldfjac, N
    REAL(SP), INTENT(IN) :: Fvec(N), X(N)
    REAL(SP), INTENT(OUT) :: Fjac(Ldfjac,N)
    !* FIRST EXECUTABLE STATEMENT  SQJAC2
    Fjac(1,1) = -1._SP
    Fjac(1,2) = 0._SP
    Fjac(2,1) = -2.E1_SP*X(1)
    Fjac(2,2) = 1.E1_SP
  END SUBROUTINE SQJAC2
  !
END MODULE TEST_SNSQE_XF
!
PROGRAM MAIN
  USE TEST_SNSQE_XF
  IMPLICIT NONE
  !
  CALL SNSQQK()
  !
END PROGRAM MAIN