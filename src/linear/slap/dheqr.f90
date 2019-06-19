!** DHEQR
SUBROUTINE DHEQR(A,Lda,N,Q,Info,Ijob)
  !> Internal routine for DGMRES.
  !***
  ! **Library:**   SLATEC (SLAP)
  !***
  ! **Category:**  D2A4, D2B4
  !***
  ! **Type:**      DOUBLE PRECISION (SHEQR-S, DHEQR-D)
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
  !        This   routine  performs  a QR   decomposition  of an  upper
  !        Hessenberg matrix A using Givens  rotations.  There  are two
  !        options  available: 1)  Performing  a fresh decomposition 2)
  !        updating the QR factors by adding a row and  a column to the
  !        matrix A.
  !
  !- Usage:
  !      INTEGER LDA, N, INFO, IJOB
  !      DOUBLE PRECISION A(LDA,N), Q(2*N)
  !
  !      CALL DHEQR(A, LDA, N, Q, INFO, IJOB)
  !
  !- Arguments:
  ! A      :INOUT    Double Precision A(LDA,N)
  !         On input, the matrix to be decomposed.
  !         On output, the upper triangular matrix R.
  !         The factorization can be written Q*A = R, where
  !         Q is a product of Givens rotations and R is upper
  !         triangular.
  ! LDA    :IN       Integer
  !         The leading dimension of the array A.
  ! N      :IN       Integer
  !         A is an (N+1) by N Hessenberg matrix.
  ! Q      :OUT      Double Precision Q(2*N)
  !         The factors c and s of each Givens rotation used
  !         in decomposing A.
  ! INFO   :OUT      Integer
  !         = 0  normal value.
  !         = K  if  A(K,K) = 0.0 .  This is not an error
  !           condition for this subroutine, but it does
  !           indicate that DHELS will divide by zero
  !           if called.
  ! IJOB   :IN       Integer
  !         = 1     means that a fresh decomposition of the
  !                 matrix A is desired.
  !         >= 2  means that the current decomposition of A
  !                 will be updated by the addition of a row
  !                 and a column.
  !
  !***
  ! **See also:**  DGMRES
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   890404  DATE WRITTEN
  !   890404  Previous REVISION DATE
  !   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
  !   890922  Numerous changes to prologue to make closer to SLATEC
  !           standard.  (FNF)
  !   890929  Numerous changes to reduce SP/DP differences.  (FNF)
  !   910411  Prologue converted to Version 4.0 format.  (BAB)
  !   910506  Made subsidiary to DGMRES.  (FNF)
  !   920511  Added complete declaration section.  (WRB)

  !     .. Scalar Arguments ..
  INTEGER :: Ijob, Info, Lda, N
  !     .. Array Arguments ..
  REAL(DP) :: A(Lda,*), Q(*)
  !     .. Local Scalars ..
  REAL(DP) :: c, s, t, t1, t2
  INTEGER :: i, iq, j, k, km1, kp1, nm1
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, SQRT
  !* FIRST EXECUTABLE STATEMENT  DHEQR
  IF( Ijob>1 ) THEN
    !   -------------------------------------------------------------------
    !         The old factorization of a will be updated.  A row and a
    !         column has been added to the matrix A.  N by N-1 is now
    !         the old size of the matrix.
    !   -------------------------------------------------------------------
    nm1 = N - 1
    !   -------------------------------------------------------------------
    !         Multiply the new column by the N previous Givens rotations.
    !   -------------------------------------------------------------------
    DO k = 1, nm1
      i = 2*(k-1) + 1
      t1 = A(k,N)
      t2 = A(k+1,N)
      c = Q(i)
      s = Q(i+1)
      A(k,N) = c*t1 - s*t2
      A(k+1,N) = s*t1 + c*t2
    END DO
    !   -------------------------------------------------------------------
    !         Complete update of decomposition by forming last Givens
    !         rotation, and multiplying it times the column
    !         vector(A(N,N),A(NP1,N)).
    !   -------------------------------------------------------------------
    Info = 0
    t1 = A(N,N)
    t2 = A(N+1,N)
    IF( t2==0.0D0 ) THEN
      c = 1
      s = 0
    ELSEIF( ABS(t2)>=ABS(t1) ) THEN
      t = t1/t2
      s = -1.0D0/SQRT(1.0D0+t*t)
      c = -s*t
    ELSE
      t = t2/t1
      c = 1.0D0/SQRT(1.0D0+t*t)
      s = -c*t
    END IF
    iq = 2*N - 1
    Q(iq) = c
    Q(iq+1) = s
    A(N,N) = c*t1 - s*t2
    IF( A(N,N)==0.0D0 ) Info = N
    RETURN
  END IF
  !   -------------------------------------------------------------------
  !         A new factorization is desired.
  !   -------------------------------------------------------------------
  !         QR decomposition without pivoting.
  !
  Info = 0
  DO k = 1, N
    km1 = k - 1
    kp1 = k + 1
    !
    !           Compute K-th column of R.
    !           First, multiply the K-th column of A by the previous
    !           K-1 Givens rotations.
    !
    IF( km1>=1 ) THEN
      DO j = 1, km1
        i = 2*(j-1) + 1
        t1 = A(j,k)
        t2 = A(j+1,k)
        c = Q(i)
        s = Q(i+1)
        A(j,k) = c*t1 - s*t2
        A(j+1,k) = s*t1 + c*t2
      END DO
    END IF
    !
    !         Compute Givens components C and S.
    !
    iq = 2*km1 + 1
    t1 = A(k,k)
    t2 = A(kp1,k)
    IF( t2==0.0D0 ) THEN
      c = 1
      s = 0
    ELSEIF( ABS(t2)>=ABS(t1) ) THEN
      t = t1/t2
      s = -1.0D0/SQRT(1.0D0+t*t)
      c = -s*t
    ELSE
      t = t2/t1
      c = 1.0D0/SQRT(1.0D0+t*t)
      s = -c*t
    END IF
    Q(iq) = c
    Q(iq+1) = s
    A(k,k) = c*t1 - s*t2
    IF( A(k,k)==0.0D0 ) Info = k
  END DO
  RETURN
  !------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
  RETURN
END SUBROUTINE DHEQR
