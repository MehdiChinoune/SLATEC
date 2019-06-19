!** CGECO
SUBROUTINE CGECO(A,Lda,N,Ipvt,Rcond,Z)
  !> Factor a matrix using Gaussian elimination and estimate
  !            the condition number of the matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2C1
  !***
  ! **Type:**      COMPLEX (SGECO-S, DGECO-D, CGECO-C)
  !***
  ! **Keywords:**  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     CGECO factors a complex matrix by Gaussian elimination
  !     and estimates the condition of the matrix.
  !
  !     If  RCOND  is not needed, CGEFA is slightly faster.
  !     To solve  A*X = B, follow CGECO By CGESL.
  !     To Compute  INVERSE(A)*C, follow CGECO by CGESL.
  !     To compute  DETERMINANT(A), follow CGECO by CGEDI.
  !     To compute  INVERSE(A), follow CGECO by CGEDI.
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA, N)
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
  !                an integer vector of pivot indices.
  !
  !        RCOND   REAL
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
  !        Z       COMPLEX(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CAXPY, CDOTC, CGEFA, CSSCAL, SCASUM

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE blas, ONLY : CAXPY, SCABS1, SCASUM

  INTEGER :: Lda, N, Ipvt(N)
  COMPLEX(SP) :: A(Lda,N), Z(N)
  REAL(SP) :: Rcond
  !
  COMPLEX(SP) :: ek, t, wk, wkm
  REAL(SP) :: anorm, s, sm, ynorm
  INTEGER :: info, j, k, kb, kp1, l
  !
  !     COMPUTE 1-NORM OF A
  !
  !* FIRST EXECUTABLE STATEMENT  CGECO
  anorm = 0.0E0
  DO j = 1, N
    anorm = MAX(anorm,SCASUM(N,A(1,j),1))
  END DO
  !
  !     FACTOR
  !
  CALL CGEFA(A,Lda,N,Ipvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
  !     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
  !     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
  !     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
  !     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  !
  !     SOLVE CTRANS(U)*W = E
  !
  ek = (1.0E0,0.0E0)
  DO j = 1, N
    Z(j) = (0.0E0,0.0E0)
  END DO
  DO k = 1, N
    IF( SCABS1(Z(k))/=0.0E0 ) ek = CSIGN1(ek,-Z(k))
    IF( SCABS1(ek-Z(k))>SCABS1(A(k,k)) ) THEN
      s = SCABS1(A(k,k))/SCABS1(ek-Z(k))
      Z = s*Z
      ek = CMPLX(s,0.0E0)*ek
    END IF
    wk = ek - Z(k)
    wkm = -ek - Z(k)
    s = SCABS1(wk)
    sm = SCABS1(wkm)
    IF( SCABS1(A(k,k))==0.0E0 ) THEN
      wk = (1.0E0,0.0E0)
      wkm = (1.0E0,0.0E0)
    ELSE
      wk = wk/CONJG(A(k,k))
      wkm = wkm/CONJG(A(k,k))
    END IF
    kp1 = k + 1
    IF( kp1<=N ) THEN
      DO j = kp1, N
        sm = sm + SCABS1(Z(j)+wkm*CONJG(A(k,j)))
        Z(j) = Z(j) + wk*CONJG(A(k,j))
        s = s + SCABS1(Z(j))
      END DO
      IF( s<sm ) THEN
        t = wkm - wk
        wk = wkm
        DO j = kp1, N
          Z(j) = Z(j) + t*CONJG(A(k,j))
        END DO
      END IF
    END IF
    Z(k) = wk
  END DO
  s = 1.0E0/SCASUM(N,Z,1)
  Z = s*Z
  !
  !     SOLVE CTRANS(L)*Y = W
  !
  DO kb = 1, N
    k = N + 1 - kb
    IF( k<N ) Z(k) = Z(k) + DOT_PRODUCT(A(k+1:N,k),Z(k+1:N))
    IF( SCABS1(Z(k))>1.0E0 ) THEN
      s = 1.0E0/SCABS1(Z(k))
      Z = s*Z
    END IF
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
  END DO
  s = 1.0E0/SCASUM(N,Z,1)
  Z = s*Z
  !
  ynorm = 1.0E0
  !
  !     SOLVE L*V = Y
  !
  DO k = 1, N
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
    IF( k<N ) CALL CAXPY(N-k,t,A(k+1,k),1,Z(k+1),1)
    IF( SCABS1(Z(k))>1.0E0 ) THEN
      s = 1.0E0/SCABS1(Z(k))
      Z = s*Z
      ynorm = s*ynorm
    END IF
  END DO
  s = 1.0E0/SCASUM(N,Z,1)
  Z = s*Z
  ynorm = s*ynorm
  !
  !     SOLVE  U*Z = V
  !
  DO kb = 1, N
    k = N + 1 - kb
    IF( SCABS1(Z(k))>SCABS1(A(k,k)) ) THEN
      s = SCABS1(A(k,k))/SCABS1(Z(k))
      Z = s*Z
      ynorm = s*ynorm
    END IF
    IF( SCABS1(A(k,k))/=0.0E0 ) Z(k) = Z(k)/A(k,k)
    IF( SCABS1(A(k,k))==0.0E0 ) Z(k) = (1.0E0,0.0E0)
    t = -Z(k)
    CALL CAXPY(k-1,t,A(1,k),1,Z(1),1)
  END DO
  !     MAKE ZNORM = 1.0
  s = 1.0E0/SCASUM(N,Z,1)
  Z = s*Z
  ynorm = s*ynorm
  !
  IF( anorm/=0.0E0 ) Rcond = ynorm/anorm
  IF( anorm==0.0E0 ) Rcond = 0.0E0
END SUBROUTINE CGECO
