!** DGECO
PURE SUBROUTINE DGECO(A,Lda,N,Ipvt,Rcond,Z)
  !> Factor a matrix using Gaussian elimination and estimate
  !            the condition number of the matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A1
  !***
  ! **Type:**      DOUBLE PRECISION (SGECO-S, DGECO-D, CGECO-C)
  !***
  ! **Keywords:**  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DGECO factors a double precision matrix by Gaussian elimination
  !     and estimates the condition of the matrix.
  !
  !     If  RCOND  is not needed, DGEFA is slightly faster.
  !     To solve  A*X = B, follow DGECO by DGESL.
  !     To compute  INVERSE(A)*C, follow DGECO by DGESL.
  !     To compute  DETERMINANT(A), follow DGECO by DGEDI.
  !     To compute  INVERSE(A), follow DGECO by DGEDI.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the matrix to be factored.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        A       an upper triangular matrix and the multipliers
  !                which were used to obtain it.
  !                The factorization can be written  A = L*U  where
  !                L  is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an INTEGER vector of pivot indices.
  !
  !        RCOND   DOUBLE PRECISION
  !                an estimate of the reciprocal condition of  A .
  !                For the system  A*X = B, relative perturbations
  !                in  A  and  B  of size  EPSILON  may cause
  !                relative perturbations in  X  of size  EPSILON/RCOND .
  !                If  RCOND  is so small that the logical expression
  !                           1.0 + RCOND = 1.0
  !                is true, then  A  may be singular to working
  !                precision.  In particular,  RCOND  is zero  if
  !                exact singularity is detected or the estimate
  !                underflows.
  !
  !        Z       DOUBLE PRECISION(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DASUM, DAXPY, DDOT, DGEFA, DSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE blas, ONLY : DAXPY

  INTEGER, INTENT(IN) :: Lda, N
  INTEGER, INTENT(OUT) :: Ipvt(N)
  REAL(DP), INTENT(INOUT) :: A(Lda,N)
  REAL(DP), INTENT(OUT) :: Rcond, Z(N)
  !
  REAL(DP) :: ek, t, wk, wkm
  REAL(DP) :: anorm, s, sm, ynorm
  INTEGER :: info, j, k, kb, kp1, l
  !
  !     COMPUTE 1-NORM OF A
  !
  !* FIRST EXECUTABLE STATEMENT  DGECO
  anorm = 0._DP
  DO j = 1, N
    anorm = MAX(anorm, SUM( ABS(A(1:N,j)) ) )
  END DO
  !
  !     FACTOR
  !
  CALL DGEFA(A,Lda,N,Ipvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
  !     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
  !     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
  !     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
  !     OVERFLOW.
  !
  !     SOLVE TRANS(U)*W = E
  !
  ek = 1._DP
  DO j = 1, N
    Z(j) = 0._DP
  END DO
  DO k = 1, N
    IF( Z(k)/=0._DP ) ek = SIGN(ek,-Z(k))
    IF( ABS(ek-Z(k))>ABS(A(k,k)) ) THEN
      s = ABS(A(k,k))/ABS(ek-Z(k))
      Z = s*Z
      ek = s*ek
    END IF
    wk = ek - Z(k)
    wkm = -ek - Z(k)
    s = ABS(wk)
    sm = ABS(wkm)
    IF( A(k,k)==0._DP ) THEN
      wk = 1._DP
      wkm = 1._DP
    ELSE
      wk = wk/A(k,k)
      wkm = wkm/A(k,k)
    END IF
    kp1 = k + 1
    IF( kp1<=N ) THEN
      DO j = kp1, N
        sm = sm + ABS(Z(j)+wkm*A(k,j))
        Z(j) = Z(j) + wk*A(k,j)
        s = s + ABS(Z(j))
      END DO
      IF( s<sm ) THEN
        t = wkm - wk
        wk = wkm
        DO j = kp1, N
          Z(j) = Z(j) + t*A(k,j)
        END DO
      END IF
    END IF
    Z(k) = wk
  END DO
  s = 1._DP/SUM( ABS(Z) )
  Z = s*Z
  !
  !     SOLVE TRANS(L)*Y = W
  !
  DO kb = 1, N
    k = N + 1 - kb
    IF( k<N ) Z(k) = Z(k) + DOT_PRODUCT(A(k+1:N,k),Z(k+1:N))
    IF( ABS(Z(k))>1._DP ) THEN
      s = 1._DP/ABS(Z(k))
      Z = s*Z
    END IF
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
  END DO
  s = 1._DP/SUM( ABS(Z) )
  Z = s*Z
  !
  ynorm = 1._DP
  !
  !     SOLVE L*V = Y
  !
  DO k = 1, N
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
    IF( k<N ) CALL DAXPY(N-k,t,A(k+1,k),1,Z(k+1),1)
    IF( ABS(Z(k))>1._DP ) THEN
      s = 1._DP/ABS(Z(k))
      Z = s*Z
      ynorm = s*ynorm
    END IF
  END DO
  s = 1._DP/SUM( ABS(Z) )
  Z = s*Z
  ynorm = s*ynorm
  !
  !     SOLVE  U*Z = V
  !
  DO kb = 1, N
    k = N + 1 - kb
    IF( ABS(Z(k))>ABS(A(k,k)) ) THEN
      s = ABS(A(k,k))/ABS(Z(k))
      Z = s*Z
      ynorm = s*ynorm
    END IF
    IF( A(k,k)/=0._DP ) Z(k) = Z(k)/A(k,k)
    IF( A(k,k)==0._DP ) Z(k) = 1._DP
    t = -Z(k)
    CALL DAXPY(k-1,t,A(1,k),1,Z(1),1)
  END DO
  !     MAKE ZNORM = 1.0
  s = 1._DP/SUM( ABS(Z) )
  Z = s*Z
  ynorm = s*ynorm
  !
  IF( anorm/=0._DP ) Rcond = ynorm/anorm
  IF( anorm==0._DP ) Rcond = 0._DP

END SUBROUTINE DGECO