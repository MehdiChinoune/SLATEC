!*==CHICO.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CHICO
SUBROUTINE CHICO(A,Lda,N,Kpvt,Rcond,Z)
  IMPLICIT NONE
  !*--CHICO5
  !***BEGIN PROLOGUE  CHICO
  !***PURPOSE  Factor a complex Hermitian matrix by elimination with sym-
  !            metric pivoting and estimate the condition of the matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D1A
  !***TYPE      COMPLEX (SSICO-S, DSICO-D, CHICO-C, CSICO-C)
  !***KEYWORDS  CONDITION NUMBER, HERMITIAN, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     CHICO factors a complex Hermitian matrix by elimination with
  !     symmetric pivoting and estimates the condition of the matrix.
  !
  !     If  RCOND  is not needed, CHIFA is slightly faster.
  !     To solve  A*X = B , follow CHICO by CHISL.
  !     To compute  INVERSE(A)*C , follow CHICO by CHISL.
  !     To compute  INVERSE(A) , follow CHICO by CHIDI.
  !     To compute  DETERMINANT(A) , follow CHICO by CHIDI.
  !     To compute  INERTIA(A), follow CHICO by CHIDI.
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA, N)
  !                the Hermitian matrix to be factored.
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
  !                The factorization can be written  A = U*D*CTRANS(U)
  !                where  U  is a product of permutation and unit
  !                upper triangular matrices , CTRANS(U) is the
  !                conjugate transpose of  U , and  D  is block diagonal
  !                with 1 by 1 and 2 by 2 blocks.
  !
  !        KVPT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        RCOND   REAL
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
  !        Z       COMPLEX(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CAXPY, CDOTC, CHIFA, CSSCAL, SCASUM
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
  !***END PROLOGUE  CHICO
  INTEGER Lda , N , Kpvt(*)
  COMPLEX A(Lda,*) , Z(*)
  REAL Rcond
  !
  COMPLEX ak , akm1 , bk , bkm1 , CDOTC , denom , ek , t
  REAL anorm , s , SCASUM , ynorm
  INTEGER i , info , j , jm1 , k , kp , kps , ks
  REAL, EXTERNAL :: CABS1
  COMPLEX, EXTERNAL :: CSIGN1
  !
  !     FIND NORM OF A USING ONLY UPPER HALF
  !
  !***FIRST EXECUTABLE STATEMENT  CHICO
  DO j = 1 , N
    Z(j) = CMPLX(SCASUM(j,A(1,j),1),0.0E0)
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1 , jm1
        Z(i) = CMPLX(REAL(Z(i))+CABS1(A(i,j)),0.0E0)
      ENDDO
    ENDIF
  ENDDO
  anorm = 0.0E0
  DO j = 1 , N
    anorm = MAX(anorm,REAL(Z(j)))
  ENDDO
  !
  !     FACTOR
  !
  CALL CHIFA(A,Lda,N,Kpvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
  !     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
  !     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
  !     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  !
  !     SOLVE U*D*W = E
  !
  ek = (1.0E0,0.0E0)
  DO j = 1 , N
    Z(j) = (0.0E0,0.0E0)
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
    IF ( CABS1(Z(k))/=0.0E0 ) ek = CSIGN1(ek,Z(k))
    Z(k) = Z(k) + ek
    CALL CAXPY(k-ks,Z(k),A(1,k),1,Z(1),1)
    IF ( ks/=1 ) THEN
      IF ( CABS1(Z(k-1))/=0.0E0 ) ek = CSIGN1(ek,Z(k-1))
      Z(k-1) = Z(k-1) + ek
      CALL CAXPY(k-ks,Z(k-1),A(1,k-1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      ak = A(k,k)/CONJG(A(k-1,k))
      akm1 = A(k-1,k-1)/A(k-1,k)
      bk = Z(k)/CONJG(A(k-1,k))
      bkm1 = Z(k-1)/A(k-1,k)
      denom = ak*akm1 - 1.0E0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( CABS1(Z(k))>CABS1(A(k,k)) ) THEN
        s = CABS1(A(k,k))/CABS1(Z(k))
        CALL CSSCAL(N,s,Z,1)
        ek = CMPLX(s,0.0E0)*ek
      ENDIF
      IF ( CABS1(A(k,k))/=0.0E0 ) Z(k) = Z(k)/A(k,k)
      IF ( CABS1(A(k,k))==0.0E0 ) Z(k) = (1.0E0,0.0E0)
    ENDIF
    k = k - ks
  ENDDO
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  !
  !     SOLVE CTRANS(U)*Y = W
  !
  k = 1
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + CDOTC(k-1,A(1,k),1,Z(1),1)
      IF ( ks==2 ) Z(k+1) = Z(k+1) + CDOTC(k-1,A(1,k+1),1,Z(1),1)
      kp = ABS(Kpvt(k))
      IF ( kp/=k ) THEN
        t = Z(k)
        Z(k) = Z(kp)
        Z(kp) = t
      ENDIF
    ENDIF
    k = k + ks
  ENDDO
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  !
  ynorm = 1.0E0
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
      CALL CAXPY(k-ks,Z(k),A(1,k),1,Z(1),1)
      IF ( ks==2 ) CALL CAXPY(k-ks,Z(k-1),A(1,k-1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      ak = A(k,k)/CONJG(A(k-1,k))
      akm1 = A(k-1,k-1)/A(k-1,k)
      bk = Z(k)/CONJG(A(k-1,k))
      bkm1 = Z(k-1)/A(k-1,k)
      denom = ak*akm1 - 1.0E0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( CABS1(Z(k))>CABS1(A(k,k)) ) THEN
        s = CABS1(A(k,k))/CABS1(Z(k))
        CALL CSSCAL(N,s,Z,1)
        ynorm = s*ynorm
      ENDIF
      IF ( CABS1(A(k,k))/=0.0E0 ) Z(k) = Z(k)/A(k,k)
      IF ( CABS1(A(k,k))==0.0E0 ) Z(k) = (1.0E0,0.0E0)
    ENDIF
    k = k - ks
  ENDDO
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  !     SOLVE CTRANS(U)*Z = V
  !
  k = 1
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + CDOTC(k-1,A(1,k),1,Z(1),1)
      IF ( ks==2 ) Z(k+1) = Z(k+1) + CDOTC(k-1,A(1,k+1),1,Z(1),1)
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
  s = 1.0E0/SCASUM(N,Z,1)
  CALL CSSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  IF ( anorm/=0.0E0 ) Rcond = ynorm/anorm
  IF ( anorm==0.0E0 ) Rcond = 0.0E0
END SUBROUTINE CHICO
