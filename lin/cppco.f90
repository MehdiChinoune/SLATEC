!DECK CPPCO
SUBROUTINE CPPCO(Ap,N,Rcond,Z,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CPPCO
  !***PURPOSE  Factor a complex Hermitian positive definite matrix stored
  !            in packed form and estimate the condition number of the
  !            matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D1B
  !***TYPE      COMPLEX (SPPCO-S, DPPCO-D, CPPCO-C)
  !***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION, PACKED, POSITIVE DEFINITE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     CPPCO factors a complex Hermitian positive definite matrix
  !     stored in packed form and estimates the condition of the matrix.
  !
  !     If  RCOND  is not needed, CPPFA is slightly faster.
  !     To solve  A*X = B, follow CPPCO by CPPSL.
  !     To compute  INVERSE(A)*C, follow CPPCO by CPPSL.
  !     To compute  DETERMINANT(A), follow CPPCO by CPPDI.
  !     To compute  INVERSE(A), follow CPPCO by CPPDI.
  !
  !     On Entry
  !
  !        AP      COMPLEX (N*(N+1)/2)
  !                the packed form of a Hermitian matrix  A .  The
  !                columns of the upper triangle are stored sequentially
  !                in a one-dimensional array of length  N*(N+1)/2 .
  !                See comments below for details.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        AP      an upper triangular matrix  R, stored in packed
  !                form, so that  A = CTRANS(R)*R .
  !                If  INFO .NE. 0, the factorization is not complete.
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
  !                underflows.  If INFO .NE. 0, RCOND is unchanged.
  !
  !        Z       COMPLEX(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is singular to working precision, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
  !                If  INFO .NE. 0, Z  is unchanged.
  !
  !        INFO    INTEGER
  !                = 0  for normal return.
  !                = K  signals an error condition.  The leading minor
  !                     of order  K  is not positive definite.
  !
  !     Packed Storage
  !
  !          The following program segment will pack the upper
  !          triangle of a Hermitian matrix.
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
  !***ROUTINES CALLED  CAXPY, CDOTC, CPPFA, CSSCAL, SCASUM
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CPPCO
  INTEGER N, Info
  COMPLEX Ap(*), Z(*)
  REAL Rcond
  !
  COMPLEX CDOTC, ek, t, wk, wkm
  REAL anorm, s, SCASUM, sm, ynorm
  INTEGER i, ij, j, jm1, j1, k, kb, kj, kk, kp1
  REAL, EXTERNAL :: CABS1
  COMPLEX, EXTERNAL :: CSIGN1
  !
  !     FIND NORM OF A
  !
  !***FIRST EXECUTABLE STATEMENT  CPPCO
  j1 = 1
  DO j = 1, N
    Z(j) = CMPLX(SCASUM(j,Ap(j1),1),0.0E0)
    ij = j1
    j1 = j1 + j
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1, jm1
        Z(i) = CMPLX(REAL(Z(i))+CABS1(Ap(ij)),0.0E0)
        ij = ij + 1
      ENDDO
    ENDIF
  ENDDO
  anorm = 0.0E0
  DO j = 1, N
    anorm = MAX(anorm,REAL(Z(j)))
  ENDDO
  !
  !     FACTOR
  !
  CALL CPPFA(Ap,N,Info)
  IF ( Info==0 ) THEN
    !
    !        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
    !        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
    !        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
    !        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
    !        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
    !
    !        SOLVE CTRANS(R)*W = E
    !
    ek = (1.0E0,0.0E0)
    DO j = 1, N
      Z(j) = (0.0E0,0.0E0)
    ENDDO
    kk = 0
    DO k = 1, N
      kk = kk + k
      IF ( CABS1(Z(k))/=0.0E0 ) ek = CSIGN1(ek,-Z(k))
      IF ( CABS1(ek-Z(k))>REAL(Ap(kk)) ) THEN
        s = REAL(Ap(kk))/CABS1(ek-Z(k))
        CALL CSSCAL(N,s,Z,1)
        ek = CMPLX(s,0.0E0)*ek
      ENDIF
      wk = ek - Z(k)
      wkm = -ek - Z(k)
      s = CABS1(wk)
      sm = CABS1(wkm)
      wk = wk/Ap(kk)
      wkm = wkm/Ap(kk)
      kp1 = k + 1
      kj = kk + k
      IF ( kp1<=N ) THEN
        DO j = kp1, N
          sm = sm + CABS1(Z(j)+wkm*CONJG(Ap(kj)))
          Z(j) = Z(j) + wk*CONJG(Ap(kj))
          s = s + CABS1(Z(j))
          kj = kj + j
        ENDDO
        IF ( s<sm ) THEN
          t = wkm - wk
          wk = wkm
          kj = kk + k
          DO j = kp1, N
            Z(j) = Z(j) + t*CONJG(Ap(kj))
            kj = kj + j
          ENDDO
        ENDIF
      ENDIF
      Z(k) = wk
    ENDDO
    s = 1.0E0/SCASUM(N,Z,1)
    CALL CSSCAL(N,s,Z,1)
    !
    !        SOLVE R*Y = W
    !
    DO kb = 1, N
      k = N + 1 - kb
      IF ( CABS1(Z(k))>REAL(Ap(kk)) ) THEN
        s = REAL(Ap(kk))/CABS1(Z(k))
        CALL CSSCAL(N,s,Z,1)
      ENDIF
      Z(k) = Z(k)/Ap(kk)
      kk = kk - k
      t = -Z(k)
      CALL CAXPY(k-1,t,Ap(kk+1),1,Z(1),1)
    ENDDO
    s = 1.0E0/SCASUM(N,Z,1)
    CALL CSSCAL(N,s,Z,1)
    !
    ynorm = 1.0E0
    !
    !        SOLVE CTRANS(R)*V = Y
    !
    DO k = 1, N
      Z(k) = Z(k) - CDOTC(k-1,Ap(kk+1),1,Z(1),1)
      kk = kk + k
      IF ( CABS1(Z(k))>REAL(Ap(kk)) ) THEN
        s = REAL(Ap(kk))/CABS1(Z(k))
        CALL CSSCAL(N,s,Z,1)
        ynorm = s*ynorm
      ENDIF
      Z(k) = Z(k)/Ap(kk)
    ENDDO
    s = 1.0E0/SCASUM(N,Z,1)
    CALL CSSCAL(N,s,Z,1)
    ynorm = s*ynorm
    !
    !        SOLVE R*Z = V
    !
    DO kb = 1, N
      k = N + 1 - kb
      IF ( CABS1(Z(k))>REAL(Ap(kk)) ) THEN
        s = REAL(Ap(kk))/CABS1(Z(k))
        CALL CSSCAL(N,s,Z,1)
        ynorm = s*ynorm
      ENDIF
      Z(k) = Z(k)/Ap(kk)
      kk = kk - k
      t = -Z(k)
      CALL CAXPY(k-1,t,Ap(kk+1),1,Z(1),1)
    ENDDO
    !        MAKE ZNORM = 1.0
    s = 1.0E0/SCASUM(N,Z,1)
    CALL CSSCAL(N,s,Z,1)
    ynorm = s*ynorm
    !
    IF ( anorm/=0.0E0 ) Rcond = ynorm/anorm
    IF ( anorm==0.0E0 ) Rcond = 0.0E0
  ENDIF
END SUBROUTINE CPPCO
