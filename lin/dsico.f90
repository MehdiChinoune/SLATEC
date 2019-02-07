!*==DSICO.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DSICO
SUBROUTINE DSICO(A,Lda,N,Kpvt,Rcond,Z)
  IMPLICIT NONE
  !*--DSICO5
  !***BEGIN PROLOGUE  DSICO
  !***PURPOSE  Factor a symmetric matrix by elimination with symmetric
  !            pivoting and estimate the condition number of the matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2B1A
  !***TYPE      DOUBLE PRECISION (SSICO-S, DSICO-D, CHICO-C, CSICO-C)
  !***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION, SYMMETRIC
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DSICO factors a double precision symmetric matrix by elimination
  !     with symmetric pivoting and estimates the condition of the
  !     matrix.
  !
  !     If  RCOND  is not needed, DSIFA is slightly faster.
  !     To solve  A*X = B , follow DSICO by DSISL.
  !     To compute  INVERSE(A)*C , follow DSICO by DSISL.
  !     To compute  INVERSE(A) , follow DSICO by DSIDI.
  !     To compute  DETERMINANT(A) , follow DSICO by DSIDI.
  !     To compute  INERTIA(A), follow DSICO by DSIDI.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the symmetric matrix to be factored.
  !                Only the diagonal and upper triangle are used.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     Output
  !
  !        A       a block diagonal matrix and the multipliers which
  !                were used to obtain it.
  !                The factorization can be written  A = U*D*TRANS(U)
  !                where  U  is a product of permutation and unit
  !                upper triangular matrices, TRANS(U) is the
  !                transpose of  U , and  D  is block diagonal
  !                with 1 by 1 and 2 by 2 blocks.
  !
  !        KPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        RCOND   DOUBLE PRECISION
  !                an estimate of the reciprocal condition of  A .
  !                For the system  A*X = B , relative perturbations
  !                in  A  and  B  of size  EPSILON  may cause
  !                relative perturbations in  X  of size  EPSILON/RCOND .
  !                If  RCOND  is so small that the logical expression
  !                           1.0 + RCOND .EQ. 1.0
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
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DASUM, DAXPY, DDOT, DSCAL, DSIFA
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Modified routine equivalence list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DSICO
  INTEGER Lda , N , Kpvt(*)
  REAL(8) :: A(Lda,*) , Z(*)
  REAL(8) :: Rcond
  !
  REAL(8) :: ak , akm1 , bk , bkm1 , DDOT , denom , ek , t
  REAL(8) :: anorm , s , DASUM , ynorm
  INTEGER i , info , j , jm1 , k , kp , kps , ks
  !
  !     FIND NORM OF A USING ONLY UPPER HALF
  !
  !***FIRST EXECUTABLE STATEMENT  DSICO
  DO j = 1 , N
    Z(j) = DASUM(j,A(1,j),1)
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1 , jm1
        Z(i) = Z(i) + ABS(A(i,j))
      ENDDO
    ENDIF
  ENDDO
  anorm = 0.0D0
  DO j = 1 , N
    anorm = MAX(anorm,Z(j))
  ENDDO
  !
  !     FACTOR
  !
  CALL DSIFA(A,Lda,N,Kpvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
  !     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
  !     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
  !     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  !
  !     SOLVE U*D*W = E
  !
  ek = 1.0D0
  DO j = 1 , N
    Z(j) = 0.0D0
  ENDDO
  k = N
  DO WHILE ( k/=0 )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    kp = ABS(Kpvt(k))
    kps = k + 1 - ks
    IF ( kp/=kps ) THEN
      t = Z(kps)
      Z(kps) = Z(kp)
      Z(kp) = t
    ENDIF
    IF ( Z(k)/=0.0D0 ) ek = SIGN(ek,Z(k))
    Z(k) = Z(k) + ek
    CALL DAXPY(k-ks,Z(k),A(1,k),1,Z(1),1)
    IF ( ks/=1 ) THEN
      IF ( Z(k-1)/=0.0D0 ) ek = SIGN(ek,Z(k-1))
      Z(k-1) = Z(k-1) + ek
      CALL DAXPY(k-ks,Z(k-1),A(1,k-1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      ak = A(k,k)/A(k-1,k)
      akm1 = A(k-1,k-1)/A(k-1,k)
      bk = Z(k)/A(k-1,k)
      bkm1 = Z(k-1)/A(k-1,k)
      denom = ak*akm1 - 1.0D0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( ABS(Z(k))>ABS(A(k,k)) ) THEN
        s = ABS(A(k,k))/ABS(Z(k))
        CALL DSCAL(N,s,Z,1)
        ek = s*ek
      ENDIF
      IF ( A(k,k)/=0.0D0 ) Z(k) = Z(k)/A(k,k)
      IF ( A(k,k)==0.0D0 ) Z(k) = 1.0D0
    ENDIF
    k = k - ks
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  !
  !     SOLVE TRANS(U)*Y = W
  !
  k = 1
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + DDOT(k-1,A(1,k),1,Z(1),1)
      IF ( ks==2 ) Z(k+1) = Z(k+1) + DDOT(k-1,A(1,k+1),1,Z(1),1)
      kp = ABS(Kpvt(k))
      IF ( kp/=k ) THEN
        t = Z(k)
        Z(k) = Z(kp)
        Z(kp) = t
      ENDIF
    ENDIF
    k = k + ks
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  !
  ynorm = 1.0D0
  !
  !     SOLVE U*D*V = Y
  !
  k = N
  DO WHILE ( k/=0 )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=ks ) THEN
      kp = ABS(Kpvt(k))
      kps = k + 1 - ks
      IF ( kp/=kps ) THEN
        t = Z(kps)
        Z(kps) = Z(kp)
        Z(kp) = t
      ENDIF
      CALL DAXPY(k-ks,Z(k),A(1,k),1,Z(1),1)
      IF ( ks==2 ) CALL DAXPY(k-ks,Z(k-1),A(1,k-1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      ak = A(k,k)/A(k-1,k)
      akm1 = A(k-1,k-1)/A(k-1,k)
      bk = Z(k)/A(k-1,k)
      bkm1 = Z(k-1)/A(k-1,k)
      denom = ak*akm1 - 1.0D0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( ABS(Z(k))>ABS(A(k,k)) ) THEN
        s = ABS(A(k,k))/ABS(Z(k))
        CALL DSCAL(N,s,Z,1)
        ynorm = s*ynorm
      ENDIF
      IF ( A(k,k)/=0.0D0 ) Z(k) = Z(k)/A(k,k)
      IF ( A(k,k)==0.0D0 ) Z(k) = 1.0D0
    ENDIF
    k = k - ks
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  !     SOLVE TRANS(U)*Z = V
  !
  k = 1
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + DDOT(k-1,A(1,k),1,Z(1),1)
      IF ( ks==2 ) Z(k+1) = Z(k+1) + DDOT(k-1,A(1,k+1),1,Z(1),1)
      kp = ABS(Kpvt(k))
      IF ( kp/=k ) THEN
        t = Z(k)
        Z(k) = Z(kp)
        Z(kp) = t
      ENDIF
    ENDIF
    k = k + ks
  ENDDO
  !     MAKE ZNORM = 1.0
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  IF ( anorm/=0.0D0 ) Rcond = ynorm/anorm
  IF ( anorm==0.0D0 ) Rcond = 0.0D0
END SUBROUTINE DSICO
