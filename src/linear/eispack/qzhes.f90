!** QZHES
SUBROUTINE QZHES(Nm,N,A,B,Matz,Z)
  !>
  !  The first step of the QZ algorithm for solving generalized
  !            matrix eigenproblems.  Accepts a pair of real general
  !            matrices and reduces one of them to upper Hessenberg
  !            and the other to upper triangular form using orthogonal
  !            transformations. Usually followed by QZIT, QZVAL, QZVEC.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C1B3
  !***
  ! **Type:**      SINGLE PRECISION (QZHES-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is the first step of the QZ algorithm
  !     for solving generalized matrix eigenvalue problems,
  !     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
  !
  !     This subroutine accepts a pair of REAL GENERAL matrices and
  !     reduces one of them to upper Hessenberg form and the other
  !     to upper triangular form using orthogonal transformations.
  !     It is usually followed by  QZIT,  QZVAL  and, possibly,  QZVEC.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A, B, and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrices A and B.  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        A contains a real general matrix.  A is a two-dimensional
  !          REAL array, dimensioned A(NM,N).
  !
  !        B contains a real general matrix.  B is a two-dimensional
  !          REAL array, dimensioned B(NM,N).
  !
  !        MATZ should be set to .TRUE. if the right hand transformations
  !          are to be accumulated for later use in computing
  !          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
  !          variable.
  !
  !     On Output
  !
  !        A has been reduced to upper Hessenberg form.  The elements
  !          below the first subdiagonal have been set to zero.
  !
  !        B has been reduced to upper triangular form.  The elements
  !          below the main diagonal have been set to zero.
  !
  !        Z contains the product of the right hand transformations if
  !          MATZ has been set to .TRUE.  Otherwise, Z is not referenced.
  !          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***
  ! **References:**  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, N, lb, l1, Nm, nk1, nm1, nm2
  REAL A(Nm,*), B(Nm,*), Z(Nm,*)
  REAL r, s, t, u1, u2, v1, v2, rho
  LOGICAL Matz
  !
  !     .......... INITIALIZE Z ..........
  !* FIRST EXECUTABLE STATEMENT  QZHES
  IF ( Matz ) THEN
    !
    DO i = 1, N
      !
      DO j = 1, N
        Z(i,j) = 0.0E0
      END DO
      !
      Z(i,i) = 1.0E0
    END DO
  END IF
  !     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
  IF ( N>1 ) THEN
    nm1 = N - 1
    !
    DO l = 1, nm1
      l1 = l + 1
      s = 0.0E0
      !
      DO i = l1, N
        s = s + ABS(B(i,l))
      END DO
      !
      IF ( s/=0.0E0 ) THEN
        s = s + ABS(B(l,l))
        r = 0.0E0
        !
        DO i = l, N
          B(i,l) = B(i,l)/s
          r = r + B(i,l)**2
        END DO
        !
        r = SIGN(SQRT(r),B(l,l))
        B(l,l) = B(l,l) + r
        rho = r*B(l,l)
        !
        DO j = l1, N
          t = 0.0E0
          !
          DO i = l, N
            t = t + B(i,l)*B(i,j)
          END DO
          !
          t = -t/rho
          !
          DO i = l, N
            B(i,j) = B(i,j) + t*B(i,l)
          END DO
          !
        END DO
        !
        DO j = 1, N
          t = 0.0E0
          !
          DO i = l, N
            t = t + B(i,l)*A(i,j)
          END DO
          !
          t = -t/rho
          !
          DO i = l, N
            A(i,j) = A(i,j) + t*B(i,l)
          END DO
          !
        END DO
        !
        B(l,l) = -s*r
        !
        DO i = l1, N
          B(i,l) = 0.0E0
        END DO
      END IF
      !
    END DO
    !     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
    !                KEEPING B TRIANGULAR ..........
    IF ( N/=2 ) THEN
      nm2 = N - 2
      !
      DO k = 1, nm2
        nk1 = nm1 - k
        !     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
        DO lb = 1, nk1
          l = N - lb
          l1 = l + 1
          !     .......... ZERO A(L+1,K) ..........
          s = ABS(A(l,k)) + ABS(A(l1,k))
          IF ( s/=0.0E0 ) THEN
            u1 = A(l,k)/s
            u2 = A(l1,k)/s
            r = SIGN(SQRT(u1*u1+u2*u2),u1)
            v1 = -(u1+r)/r
            v2 = -u2/r
            u2 = v2/v1
            !
            DO j = k, N
              t = A(l,j) + u2*A(l1,j)
              A(l,j) = A(l,j) + t*v1
              A(l1,j) = A(l1,j) + t*v2
            END DO
            !
            A(l1,k) = 0.0E0
            !
            DO j = l, N
              t = B(l,j) + u2*B(l1,j)
              B(l,j) = B(l,j) + t*v1
              B(l1,j) = B(l1,j) + t*v2
            END DO
            !     .......... ZERO B(L+1,L) ..........
            s = ABS(B(l1,l1)) + ABS(B(l1,l))
            IF ( s/=0.0E0 ) THEN
              u1 = B(l1,l1)/s
              u2 = B(l1,l)/s
              r = SIGN(SQRT(u1*u1+u2*u2),u1)
              v1 = -(u1+r)/r
              v2 = -u2/r
              u2 = v2/v1
              !
              DO i = 1, l1
                t = B(i,l1) + u2*B(i,l)
                B(i,l1) = B(i,l1) + t*v1
                B(i,l) = B(i,l) + t*v2
              END DO
              !
              B(l1,l) = 0.0E0
              !
              DO i = 1, N
                t = A(i,l1) + u2*A(i,l)
                A(i,l1) = A(i,l1) + t*v1
                A(i,l) = A(i,l) + t*v2
              END DO
              !
              IF ( Matz ) THEN
                !
                DO i = 1, N
                  t = Z(i,l1) + u2*Z(i,l)
                  Z(i,l1) = Z(i,l1) + t*v1
                  Z(i,l) = Z(i,l) + t*v2
                END DO
              END IF
            END IF
          END IF
          !
        END DO
        !
      END DO
    END IF
  END IF
  !
END SUBROUTINE QZHES
