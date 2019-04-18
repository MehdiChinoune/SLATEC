!** MINFIT
SUBROUTINE MINFIT(Nm,M,N,A,W,Ip,B,Ierr,Rv1)
  !>
  !  Compute the singular value decomposition of a rectangular
  !            matrix and solve the related linear least squares problem.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D9
  !***
  ! **Type:**      SINGLE PRECISION (MINFIT-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure MINFIT,
  !     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
  !     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
  !
  !     This subroutine determines, towards the solution of the linear
  !                                                        T
  !     system AX=B, the singular value decomposition A=USV  of a real
  !                                         T
  !     M by N rectangular matrix, forming U B rather than U.  Householder
  !     bidiagonalization and a variant of the QR algorithm are used.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and B, as declared in the calling
  !          program dimension statement.  Note that NM must be at least
  !          as large as the maximum of M and N.  NM is an INTEGER
  !          variable.
  !
  !        M is the number of rows of A and B.  M is an INTEGER variable.
  !
  !        N is the number of columns of A and the order of V.  N is an
  !          INTEGER variable.
  !
  !        A contains the rectangular coefficient matrix of the system.
  !          A is a two-dimensional REAL array, dimensioned A(NM,N).
  !
  !        IP is the number of columns of B.  IP can be zero.
  !
  !        B contains the constant column matrix of the system if IP is
  !          not zero.  Otherwise, B is not referenced.  B is a two-
  !          dimensional REAL array, dimensioned B(NM,IP).
  !
  !     On OUTPUT
  !
  !        A has been overwritten by the matrix V (orthogonal) of the
  !          decomposition in its first N rows and columns.  If an
  !          error exit is made, the columns of V corresponding to
  !          indices of correct singular values should be correct.
  !
  !        W contains the N (non-negative) singular values of A (the
  !          diagonal elements of S).  They are unordered.  If an
  !          error exit is made, the singular values should be correct
  !          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
  !          REAL array, dimensioned W(N).
  !
  !                                   T
  !        B has been overwritten by U B.  If an error exit is made,
  !                       T
  !          the rows of U B corresponding to indices of correct singular
  !          values should be correct.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          K          if the K-th singular value has not been
  !                     determined after 30 iterations.
  !                     The singular values should be correct for
  !                     indices IERR+1, IERR+2, ..., N.
  !
  !        RV1 is a one-dimensional REAL array used for temporary storage,
  !          dimensioned RV1(N).
  !
  !     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
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
  ! **Routines called:**  PYTHAG

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, M, N, ii, Ip, i1, kk, k1, ll, l1, m1, Nm, its, Ierr
  REAL A(Nm,*), W(*), B(Nm,Ip), Rv1(*)
  REAL c, f, g, h, s, x, y, z, scalee, s1
  !
  !* FIRST EXECUTABLE STATEMENT  MINFIT
  Ierr = 0
  !     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
  g = 0.0E0
  scalee = 0.0E0
  s1 = 0.0E0
  !
  DO i = 1, N
    l = i + 1
    Rv1(i) = scalee*g
    g = 0.0E0
    s = 0.0E0
    scalee = 0.0E0
    IF ( i<=M ) THEN
      !
      DO k = i, M
        scalee = scalee + ABS(A(k,i))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        !
        DO k = i, M
          A(k,i) = A(k,i)/scalee
          s = s + A(k,i)**2
        END DO
        !
        f = A(i,i)
        g = -SIGN(SQRT(s),f)
        h = f*g - s
        A(i,i) = f - g
        IF ( i/=N ) THEN
          !
          DO j = l, N
            s = 0.0E0
            !
            DO k = i, M
              s = s + A(k,i)*A(k,j)
            END DO
            !
            f = s/h
            !
            DO k = i, M
              A(k,j) = A(k,j) + f*A(k,i)
            END DO
          END DO
        END IF
        !
        IF ( Ip/=0 ) THEN
          !
          DO j = 1, Ip
            s = 0.0E0
            !
            DO k = i, M
              s = s + A(k,i)*B(k,j)
            END DO
            !
            f = s/h
            !
            DO k = i, M
              B(k,j) = B(k,j) + f*A(k,i)
            END DO
          END DO
        END IF
        !
        DO k = i, M
          A(k,i) = scalee*A(k,i)
        END DO
      END IF
    END IF
    !
    W(i) = scalee*g
    g = 0.0E0
    s = 0.0E0
    scalee = 0.0E0
    IF ( i<=M.AND.i/=N ) THEN
      !
      DO k = l, N
        scalee = scalee + ABS(A(i,k))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        !
        DO k = l, N
          A(i,k) = A(i,k)/scalee
          s = s + A(i,k)**2
        END DO
        !
        f = A(i,l)
        g = -SIGN(SQRT(s),f)
        h = f*g - s
        A(i,l) = f - g
        !
        DO k = l, N
          Rv1(k) = A(i,k)/h
        END DO
        !
        IF ( i/=M ) THEN
          !
          DO j = l, M
            s = 0.0E0
            !
            DO k = l, N
              s = s + A(j,k)*A(i,k)
            END DO
            !
            DO k = l, N
              A(j,k) = A(j,k) + s*Rv1(k)
            END DO
          END DO
        END IF
        !
        DO k = l, N
          A(i,k) = scalee*A(i,k)
        END DO
      END IF
    END IF
    !
    s1 = MAX(s1,ABS(W(i))+ABS(Rv1(i)))
  END DO
  !     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
  !                FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO ii = 1, N
    i = N + 1 - ii
    IF ( i/=N ) THEN
      IF ( g/=0.0E0 ) THEN
        !
        DO j = l, N
          !     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
          A(j,i) = (A(i,j)/A(i,l))/g
        END DO
        !
        DO j = l, N
          s = 0.0E0
          !
          DO k = l, N
            s = s + A(i,k)*A(k,j)
          END DO
          !
          DO k = l, N
            A(k,j) = A(k,j) + s*A(k,i)
          END DO
        END DO
      END IF
      !
      DO j = l, N
        A(i,j) = 0.0E0
        A(j,i) = 0.0E0
      END DO
    END IF
    !
    A(i,i) = 1.0E0
    g = Rv1(i)
    l = i
  END DO
  !
  IF ( M<N.AND.Ip/=0 ) THEN
    m1 = M + 1
    !
    DO i = m1, N
      !
      DO j = 1, Ip
        B(i,j) = 0.0E0
      END DO
    END DO
  END IF
  !     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  !     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
  DO kk = 1, N
    k1 = N - kk
    k = k1 + 1
    its = 0
    !     .......... TEST FOR SPLITTING.
    !                FOR L=K STEP -1 UNTIL 1 DO -- ..........
    50 CONTINUE
    DO ll = 1, k
      l1 = k - ll
      l = l1 + 1
      IF ( s1+ABS(Rv1(l))==s1 ) GOTO 100
      !     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
      !                THROUGH THE BOTTOM OF THE LOOP ..........
      IF ( s1+ABS(W(l1))==s1 ) EXIT
    END DO
    !     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
    c = 0.0E0
    s = 1.0E0
    !
    DO i = l, k
      f = s*Rv1(i)
      Rv1(i) = c*Rv1(i)
      IF ( s1+ABS(f)==s1 ) EXIT
      g = W(i)
      h = PYTHAG(f,g)
      W(i) = h
      c = g/h
      s = -f/h
      IF ( Ip/=0 ) THEN
        !
        DO j = 1, Ip
          y = B(l1,j)
          z = B(i,j)
          B(l1,j) = y*c + z*s
          B(i,j) = -y*s + z*c
        END DO
      END IF
      !
    END DO
    !     .......... TEST FOR CONVERGENCE ..........
    100  z = W(k)
    IF ( l/=k ) THEN
      !     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
      IF ( its==30 ) GOTO 200
      its = its + 1
      x = W(l)
      y = W(k1)
      g = Rv1(k1)
      h = Rv1(k)
      f = 0.5E0*(((g+z)/h)*((g-z)/y)+y/h-h/y)
      g = PYTHAG(f,1.0E0)
      f = x - (z/x)*z + (h/x)*(y/(f+SIGN(g,f))-h)
      !     .......... NEXT QR TRANSFORMATION ..........
      c = 1.0E0
      s = 1.0E0
      !
      DO i1 = l, k1
        i = i1 + 1
        g = Rv1(i)
        y = W(i)
        h = s*g
        g = c*g
        z = PYTHAG(f,h)
        Rv1(i1) = z
        c = f/z
        s = h/z
        f = x*c + g*s
        g = -x*s + g*c
        h = y*s
        y = y*c
        !
        DO j = 1, N
          x = A(j,i1)
          z = A(j,i)
          A(j,i1) = x*c + z*s
          A(j,i) = -x*s + z*c
        END DO
        !
        z = PYTHAG(f,h)
        W(i1) = z
        !     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
        IF ( z/=0.0E0 ) THEN
          c = f/z
          s = h/z
        END IF
        f = c*g + s*y
        x = -s*g + c*y
        IF ( Ip/=0 ) THEN
          !
          DO j = 1, Ip
            y = B(i1,j)
            z = B(i,j)
            B(i1,j) = y*c + z*s
            B(i,j) = -y*s + z*c
          END DO
        END IF
        !
      END DO
      !
      Rv1(l) = 0.0E0
      Rv1(k) = f
      W(k) = x
      GOTO 50
      !     .......... CONVERGENCE ..........
    ELSEIF ( z<0.0E0 ) THEN
      !     .......... W(K) IS MADE NON-NEGATIVE ..........
      W(k) = -z
      !
      DO j = 1, N
        A(j,k) = -A(j,k)
      END DO
    END IF
    !
  END DO
  !
  RETURN
  !     .......... SET ERROR -- NO CONVERGENCE TO A
  !                SINGULAR VALUE AFTER 30 ITERATIONS ..........
  200  Ierr = k
  RETURN
END SUBROUTINE MINFIT
