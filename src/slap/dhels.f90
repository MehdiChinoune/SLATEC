!** DHELS
SUBROUTINE DHELS(A,Lda,N,Q,B)
  IMPLICIT NONE
  !>
  !***
  !  Internal routine for DGMRES.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (SHELS-S, DHELS-D)
  !***
  ! **Keywords:**  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
  !             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
  !***
  ! **Author:**  Brown, Peter, (LLNL), pnbrown@llnl.gov
  !           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
  !           Seager, Mark K., (LLNL), seager@llnl.gov
  !             Lawrence Livermore National Laboratory
  !             PO Box 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !***
  ! **Description:**
  !        This routine is extracted from the LINPACK routine SGESL with
  !        changes due to the fact that A is an upper Hessenberg matrix.
  !
  !        DHELS solves the least squares problem:
  !
  !                   MIN(B-A*X,B-A*X)
  !
  !        using the factors computed by DHEQR.
  !
  !- Usage:
  !      INTEGER LDA, N
  !      DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)
  !
  !      CALL DHELS(A, LDA, N, Q, B)
  !
  !- Arguments:
  ! A       :IN       Double Precision A(LDA,N)
  !          The output from DHEQR which contains the upper
  !          triangular factor R in the QR decomposition of A.
  ! LDA     :IN       Integer
  !          The leading dimension of the array A.
  ! N       :IN       Integer
  !          A is originally an (N+1) by N matrix.
  ! Q       :IN       Double Precision Q(2*N)
  !          The coefficients of the N Givens rotations
  !          used in the QR factorization of A.
  ! B       :INOUT    Double Precision B(N+1)
  !          On input, B is the right hand side vector.
  !          On output, B is the solution vector X.
  !
  !***
  ! **See also:**  DGMRES
  !***
  ! **Routines called:**  DAXPY

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
  !   910506  Made subsidiary to DGMRES.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  
  !         The following is for optimized compilation on LLNL/LTSS Crays.
  !LLL. OPTIMIZE
  !     .. Scalar Arguments ..
  INTEGER Lda, N
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), B(*), Q(*)
  !     .. Local Scalars ..
  REAL(8) :: c, s, t, t1, t2
  INTEGER iq, k, kb, kp1
  !     .. External Subroutines ..
  EXTERNAL :: DAXPY
  !* FIRST EXECUTABLE STATEMENT  DHELS
  !
  !         Minimize(B-A*X,B-A*X).  First form Q*B.
  !
  DO k = 1, N
    kp1 = k + 1
    iq = 2*(k-1) + 1
    c = Q(iq)
    s = Q(iq+1)
    t1 = B(k)
    t2 = B(kp1)
    B(k) = c*t1 - s*t2
    B(kp1) = s*t1 + c*t2
  END DO
  !
  !         Now solve  R*X = Q*B.
  !
  DO kb = 1, N
    k = N + 1 - kb
    B(k) = B(k)/A(k,k)
    t = -B(k)
    CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
  END DO
  !------------- LAST LINE OF DHELS FOLLOWS ----------------------------
END SUBROUTINE DHELS
