!** DPPCO
SUBROUTINE DPPCO(Ap,N,Rcond,Z,Info)
  !>
  !  Factor a symmetric positive definite matrix stored in
  !            packed form and estimate the condition number of the
  !            matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B
  !***
  ! **Type:**      DOUBLE PRECISION (SPPCO-S, DPPCO-D, CPPCO-C)
  !***
  ! **Keywords:**  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION, PACKED, POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPPCO factors a double precision symmetric positive definite
  !     matrix stored in packed form
  !     and estimates the condition of the matrix.
  !
  !     If  RCOND  is not needed, DPPFA is slightly faster.
  !     To solve  A*X = B, follow DPPCO by DPPSL.
  !     To compute  INVERSE(A)*C, follow DPPCO by DPPSL.
  !     To compute  DETERMINANT(A), follow DPPCO by DPPDI.
  !     To compute  INVERSE(A), follow DPPCO by DPPDI.
  !
  !     On Entry
  !
  !        AP      DOUBLE PRECISION (N*(N+1)/2)
  !                the packed form of a symmetric matrix  A .  The
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
  !                form, so that  A = TRANS(R)*R .
  !                If  INFO .NE. 0, the factorization is not complete.
  !
  !        RCOND   DOUBLE PRECISION
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
  !        Z       DOUBLE PRECISION(N)
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
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DASUM, DAXPY, DDOT, DPPFA, DSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N, Info
  REAL(DP) :: Ap(*), Z(*)
  REAL(DP) :: Rcond
  !
  REAL(DP) :: ek, t, wk, wkm
  REAL(DP) :: anorm, s, sm, ynorm
  INTEGER i, ij, j, jm1, j1, k, kb, kj, kk, kp1
  !
  !     FIND NORM OF A
  !
  !* FIRST EXECUTABLE STATEMENT  DPPCO
  j1 = 1
  DO j = 1, N
    Z(j) = DASUM(j,Ap(j1),1)
    ij = j1
    j1 = j1 + j
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1, jm1
        Z(i) = Z(i) + ABS(Ap(ij))
        ij = ij + 1
      END DO
    END IF
  END DO
  anorm = 0.0D0
  DO j = 1, N
    anorm = MAX(anorm,Z(j))
  END DO
  !
  !     FACTOR
  !
  CALL DPPFA(Ap,N,Info)
  IF ( Info==0 ) THEN
    !
    !        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
    !        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
    !        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
    !        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
    !        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
    !
    !        SOLVE TRANS(R)*W = E
    !
    ek = 1.0D0
    DO j = 1, N
      Z(j) = 0.0D0
    END DO
    kk = 0
    DO k = 1, N
      kk = kk + k
      IF ( Z(k)/=0.0D0 ) ek = SIGN(ek,-Z(k))
      IF ( ABS(ek-Z(k))>Ap(kk) ) THEN
        s = Ap(kk)/ABS(ek-Z(k))
        CALL DSCAL(N,s,Z,1)
        ek = s*ek
      END IF
      wk = ek - Z(k)
      wkm = -ek - Z(k)
      s = ABS(wk)
      sm = ABS(wkm)
      wk = wk/Ap(kk)
      wkm = wkm/Ap(kk)
      kp1 = k + 1
      kj = kk + k
      IF ( kp1<=N ) THEN
        DO j = kp1, N
          sm = sm + ABS(Z(j)+wkm*Ap(kj))
          Z(j) = Z(j) + wk*Ap(kj)
          s = s + ABS(Z(j))
          kj = kj + j
        END DO
        IF ( s<sm ) THEN
          t = wkm - wk
          wk = wkm
          kj = kk + k
          DO j = kp1, N
            Z(j) = Z(j) + t*Ap(kj)
            kj = kj + j
          END DO
        END IF
      END IF
      Z(k) = wk
    END DO
    s = 1.0D0/DASUM(N,Z,1)
    CALL DSCAL(N,s,Z,1)
    !
    !        SOLVE R*Y = W
    !
    DO kb = 1, N
      k = N + 1 - kb
      IF ( ABS(Z(k))>Ap(kk) ) THEN
        s = Ap(kk)/ABS(Z(k))
        CALL DSCAL(N,s,Z,1)
      END IF
      Z(k) = Z(k)/Ap(kk)
      kk = kk - k
      t = -Z(k)
      CALL DAXPY(k-1,t,Ap(kk+1),1,Z(1),1)
    END DO
    s = 1.0D0/DASUM(N,Z,1)
    CALL DSCAL(N,s,Z,1)
    !
    ynorm = 1.0D0
    !
    !        SOLVE TRANS(R)*V = Y
    !
    DO k = 1, N
      Z(k) = Z(k) - DDOT(k-1,Ap(kk+1),1,Z(1),1)
      kk = kk + k
      IF ( ABS(Z(k))>Ap(kk) ) THEN
        s = Ap(kk)/ABS(Z(k))
        CALL DSCAL(N,s,Z,1)
        ynorm = s*ynorm
      END IF
      Z(k) = Z(k)/Ap(kk)
    END DO
    s = 1.0D0/DASUM(N,Z,1)
    CALL DSCAL(N,s,Z,1)
    ynorm = s*ynorm
    !
    !        SOLVE R*Z = V
    !
    DO kb = 1, N
      k = N + 1 - kb
      IF ( ABS(Z(k))>Ap(kk) ) THEN
        s = Ap(kk)/ABS(Z(k))
        CALL DSCAL(N,s,Z,1)
        ynorm = s*ynorm
      END IF
      Z(k) = Z(k)/Ap(kk)
      kk = kk - k
      t = -Z(k)
      CALL DAXPY(k-1,t,Ap(kk+1),1,Z(1),1)
    END DO
    !        MAKE ZNORM = 1.0
    s = 1.0D0/DASUM(N,Z,1)
    CALL DSCAL(N,s,Z,1)
    ynorm = s*ynorm
    !
    IF ( anorm/=0.0D0 ) Rcond = ynorm/anorm
    IF ( anorm==0.0D0 ) Rcond = 0.0D0
  END IF
END SUBROUTINE DPPCO
