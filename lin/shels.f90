!*==SHELS.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SHELS
SUBROUTINE SHELS(A,Lda,N,Q,B)
  IMPLICIT NONE
  !*--SHELS5
  !***BEGIN PROLOGUE  SHELS
  !***SUBSIDIARY
  !***PURPOSE  Internal routine for SGMRES.
  !***LIBRARY   SLATEC (SLAP)
  !***CATEGORY  D2A4, D2B4
  !***TYPE      SINGLE PRECISION (SHELS-S, DHELS-D)
  !***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
  !             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
  !***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
  !           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
  !           Seager, Mark K., (LLNL), seager@llnl.gov
  !             Lawrence Livermore National Laboratory
  !             PO Box 808, L-60
  !             Livermore, CA 94550 (510) 423-3141
  !***DESCRIPTION
  !        This routine is extracted from the LINPACK routine SGESL with
  !        changes due to the fact that A is an upper Hessenberg matrix.
  !
  !        SHELS solves the least squares problem:
  !
  !                   MIN(B-A*X,B-A*X)
  !
  !        using the factors computed by SHEQR.
  !
  ! *Usage:
  !      INTEGER LDA, N
  !      REAL A(LDA,N), Q(2*N), B(N+1)
  !
  !      CALL SHELS(A, LDA, N, Q, B)
  !
  ! *Arguments:
  ! A       :IN       Real A(LDA,N)
  !          The output from SHEQR which contains the upper
  !          triangular factor R in the QR decomposition of A.
  ! LDA     :IN       Integer
  !          The leading dimension of the array A.
  ! N       :IN       Integer
  !          A is originally an (N+1) by N matrix.
  ! Q       :IN       Real Q(2*N)
  !          The coefficients of the N Givens rotations
  !          used in the QR factorization of A.
  ! B       :INOUT    Real B(N+1)
  !          On input, B is the right hand side vector.
  !          On output, B is the solution vector X.
  !
  !***SEE ALSO  SGMRES
  !***ROUTINES CALLED  SAXPY
  !***REVISION HISTORY  (YYMMDD)
  !   871001  DATE WRITTEN
  !   881213  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
  !   910506  Made subsidiary to SGMRES.  (FNF)
  !   920511  Added complete declaration section.  (WRB)
  !***END PROLOGUE  SHELS
  !         The following is for optimized compilation on LLNL/LTSS Crays.
  !LLL. OPTIMIZE
  !     .. Scalar Arguments ..
  INTEGER Lda, N
  !     .. Array Arguments ..
  REAL A(Lda,*), B(*), Q(*)
  !     .. Local Scalars ..
  REAL c, s, t, t1, t2
  INTEGER iq, k, kb, kp1
  !     .. External Subroutines ..
  EXTERNAL SAXPY
  !***FIRST EXECUTABLE STATEMENT  SHELS
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
  ENDDO
  !
  !         Now solve  R*X = Q*B.
  !
  DO kb = 1, N
    k = N + 1 - kb
    B(k) = B(k)/A(k,k)
    t = -B(k)
    CALL SAXPY(k-1,t,A(1,k),1,B(1),1)
  ENDDO
  !------------- LAST LINE OF SHELS FOLLOWS ----------------------------
END SUBROUTINE SHELS
