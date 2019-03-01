!DECK SSPCO
SUBROUTINE SSPCO(Ap,N,Kpvt,Rcond,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SSPCO
  !***PURPOSE  Factor a real symmetric matrix stored in packed form
  !            by elimination with symmetric pivoting and estimate the
  !            condition number of the matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2B1A
  !***TYPE      SINGLE PRECISION (SSPCO-S, DSPCO-D, CHPCO-C, CSPCO-C)
  !***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION, PACKED, SYMMETRIC
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     SSPCO factors a real symmetric matrix stored in packed
  !     form by elimination with symmetric pivoting and estimates
  !     the condition of the matrix.
  !
  !     If  RCOND  is not needed, SSPFA is slightly faster.
  !     To solve  A*X = B, follow SSPCO by SSPSL.
  !     To compute  INVERSE(A)*C, follow SSPCO by SSPSL.
  !     To compute  INVERSE(A), follow SSPCO by SSPDI.
  !     To compute  DETERMINANT(A), follow SSPCO by SSPDI.
  !     To compute  INERTIA(A), follow SSPCO by SSPDI.
  !
  !     On Entry
  !
  !        AP      REAL (N*(N+1)/2)
  !                the packed form of a symmetric matrix  A .  The
  !                columns of the upper triangle are stored sequentially
  !                in a one-dimensional array of length  N*(N+1)/2 .
  !                See comments below for details.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     Output
  !
  !        AP      a block diagonal matrix and the multipliers which
  !                were used to obtain it stored in packed form.
  !                The factorization can be written  A = U*D*TRANS(U)
  !                where  U  is a product of permutation and unit
  !                upper triangular matrices, TRANS(U) is the
  !                transpose of  U, and  D  is block diagonal
  !                with 1 by 1 and 2 by 2 blocks.
  !
  !        KPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        RCOND   REAL
  !                an estimate of the reciprocal condition of  A .
  !                For the system  A*X = B, relative perturbations
  !                in  A  and  B  of size  EPSILON  may cause
  !                relative perturbations in  X  of size  EPSILON/RCOND .
  !                If  RCOND  is so small that the logical expression
  !                           1.0 + RCOND .EQ. 1.0
  !                is true, then  A  may be singular to working
  !                precision.  In particular,  RCOND  is zero  if
  !                exact singularity is detected or the estimate
  !                underflows.
  !
  !        Z       REAL(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !
  !     Packed Storage
  !
  !          The following program segment will pack the upper
  !          triangle of a symmetric matrix.
  !
  !                K = 0
  !                DO 20 J = 1, N
  !                   DO 10 I = 1, J
  !                      K = K + 1
  !                      AP(K) = A(I,J)
  !             10    CONTINUE
  !             20 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  SASUM, SAXPY, SDOT, SSCAL, SSPFA
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
  !***END PROLOGUE  SSPCO
  INTEGER N, Kpvt(*)
  REAL Ap(*), Z(*)
  REAL Rcond
  !
  REAL ak, akm1, bk, bkm1, SDOT, denom, ek, t
  REAL anorm, s, SASUM, ynorm
  INTEGER i, ij, ik, ikm1, ikp1, info, j, jm1, j1
  INTEGER k, kk, km1k, km1km1, kp, kps, ks
  !
  !     FIND NORM OF A USING ONLY UPPER HALF
  !
  !***FIRST EXECUTABLE STATEMENT  SSPCO
  j1 = 1
  DO j = 1, N
    Z(j) = SASUM(j,Ap(j1),1)
    ij = j1
    j1 = j1 + j
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1, jm1
        Z(i) = Z(i) + ABS(Ap(ij))
        ij = ij + 1
      ENDDO
    ENDIF
  ENDDO
  anorm = 0.0E0
  DO j = 1, N
    anorm = MAX(anorm,Z(j))
  ENDDO
  !
  !     FACTOR
  !
  CALL SSPFA(Ap,N,Kpvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
  !     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
  !     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
  !     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
  !
  !     SOLVE U*D*W = E
  !
  ek = 1.0E0
  DO j = 1, N
    Z(j) = 0.0E0
  ENDDO
  k = N
  ik = (N*(N-1))/2
  DO WHILE ( k/=0 )
    kk = ik + k
    ikm1 = ik - (k-1)
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    kp = ABS(Kpvt(k))
    kps = k + 1 - ks
    IF ( kp/=kps ) THEN
      t = Z(kps)
      Z(kps) = Z(kp)
      Z(kp) = t
    ENDIF
    IF ( Z(k)/=0.0E0 ) ek = SIGN(ek,Z(k))
    Z(k) = Z(k) + ek
    CALL SAXPY(k-ks,Z(k),Ap(ik+1),1,Z(1),1)
    IF ( ks/=1 ) THEN
      IF ( Z(k-1)/=0.0E0 ) ek = SIGN(ek,Z(k-1))
      Z(k-1) = Z(k-1) + ek
      CALL SAXPY(k-ks,Z(k-1),Ap(ikm1+1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = Ap(kk)/Ap(km1k)
      akm1 = Ap(km1km1)/Ap(km1k)
      bk = Z(k)/Ap(km1k)
      bkm1 = Z(k-1)/Ap(km1k)
      denom = ak*akm1 - 1.0E0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( ABS(Z(k))>ABS(Ap(kk)) ) THEN
        s = ABS(Ap(kk))/ABS(Z(k))
        CALL SSCAL(N,s,Z,1)
        ek = s*ek
      ENDIF
      IF ( Ap(kk)/=0.0E0 ) Z(k) = Z(k)/Ap(kk)
      IF ( Ap(kk)==0.0E0 ) Z(k) = 1.0E0
    ENDIF
    k = k - ks
    ik = ik - k
    IF ( ks==2 ) ik = ik - (k+1)
  ENDDO
  s = 1.0E0/SASUM(N,Z,1)
  CALL SSCAL(N,s,Z,1)
  !
  !     SOLVE TRANS(U)*Y = W
  !
  k = 1
  ik = 0
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + SDOT(k-1,Ap(ik+1),1,Z(1),1)
      ikp1 = ik + k
      IF ( ks==2 ) Z(k+1) = Z(k+1) + SDOT(k-1,Ap(ikp1+1),1,Z(1),1)
      kp = ABS(Kpvt(k))
      IF ( kp/=k ) THEN
        t = Z(k)
        Z(k) = Z(kp)
        Z(kp) = t
      ENDIF
    ENDIF
    ik = ik + k
    IF ( ks==2 ) ik = ik + (k+1)
    k = k + ks
  ENDDO
  s = 1.0E0/SASUM(N,Z,1)
  CALL SSCAL(N,s,Z,1)
  !
  ynorm = 1.0E0
  !
  !     SOLVE U*D*V = Y
  !
  k = N
  ik = N*(N-1)/2
  DO WHILE ( k/=0 )
    kk = ik + k
    ikm1 = ik - (k-1)
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
      CALL SAXPY(k-ks,Z(k),Ap(ik+1),1,Z(1),1)
      IF ( ks==2 ) CALL SAXPY(k-ks,Z(k-1),Ap(ikm1+1),1,Z(1),1)
    ENDIF
    IF ( ks==2 ) THEN
      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = Ap(kk)/Ap(km1k)
      akm1 = Ap(km1km1)/Ap(km1k)
      bk = Z(k)/Ap(km1k)
      bkm1 = Z(k-1)/Ap(km1k)
      denom = ak*akm1 - 1.0E0
      Z(k) = (akm1*bk-bkm1)/denom
      Z(k-1) = (ak*bkm1-bk)/denom
    ELSE
      IF ( ABS(Z(k))>ABS(Ap(kk)) ) THEN
        s = ABS(Ap(kk))/ABS(Z(k))
        CALL SSCAL(N,s,Z,1)
        ynorm = s*ynorm
      ENDIF
      IF ( Ap(kk)/=0.0E0 ) Z(k) = Z(k)/Ap(kk)
      IF ( Ap(kk)==0.0E0 ) Z(k) = 1.0E0
    ENDIF
    k = k - ks
    ik = ik - k
    IF ( ks==2 ) ik = ik - (k+1)
  ENDDO
  s = 1.0E0/SASUM(N,Z,1)
  CALL SSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  !     SOLVE TRANS(U)*Z = V
  !
  k = 1
  ik = 0
  DO WHILE ( k<=N )
    ks = 1
    IF ( Kpvt(k)<0 ) ks = 2
    IF ( k/=1 ) THEN
      Z(k) = Z(k) + SDOT(k-1,Ap(ik+1),1,Z(1),1)
      ikp1 = ik + k
      IF ( ks==2 ) Z(k+1) = Z(k+1) + SDOT(k-1,Ap(ikp1+1),1,Z(1),1)
      kp = ABS(Kpvt(k))
      IF ( kp/=k ) THEN
        t = Z(k)
        Z(k) = Z(kp)
        Z(kp) = t
      ENDIF
    ENDIF
    ik = ik + k
    IF ( ks==2 ) ik = ik + (k+1)
    k = k + ks
  ENDDO
  !     MAKE ZNORM = 1.0
  s = 1.0E0/SASUM(N,Z,1)
  CALL SSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  IF ( anorm/=0.0E0 ) Rcond = ynorm/anorm
  IF ( anorm==0.0E0 ) Rcond = 0.0E0
END SUBROUTINE SSPCO
