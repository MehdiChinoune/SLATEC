!** SVD
SUBROUTINE SVD(Nm,M,N,A,W,Matu,U,Matv,V,Ierr,Rv1)
  !>
  !***
  !  Perform the singular value decomposition of a rectangular
  !            matrix.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SVD-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure SVD,
  !     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
  !     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
  !
  !     This subroutine determines the singular value decomposition
  !          T
  !     A=USV  of a REAL M by N rectangular matrix.  Householder
  !     bidiagonalization and a variant of the QR algorithm are used.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A, U and V, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !          Note that NM must be at least as large as the maximum
  !          of M and N.
  !
  !        M is the number of rows of A and U.
  !
  !        N is the number of columns of A and U and the order of V.
  !
  !        A contains the rectangular input matrix to be decomposed.  A is
  !          a two-dimensional REAL array, dimensioned A(NM,N).
  !
  !        MATU should be set to .TRUE. if the U matrix in the
  !          decomposition is desired, and to .FALSE. otherwise.
  !          MATU is a LOGICAL variable.
  !
  !        MATV should be set to .TRUE. if the V matrix in the
  !          decomposition is desired, and to .FALSE. otherwise.
  !          MATV is a LOGICAL variable.
  !
  !     On Output
  !
  !        A is unaltered (unless overwritten by U or V).
  !
  !        W contains the N (non-negative) singular values of A (the
  !          diagonal elements of S).  They are unordered.  If an
  !          error exit is made, the singular values should be correct
  !          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
  !          REAL array, dimensioned W(N).
  !
  !        U contains the matrix U (orthogonal column vectors) of the
  !          decomposition if MATU has been set to .TRUE.  Otherwise,
  !          U is used as a temporary array.  U may coincide with A.
  !          If an error exit is made, the columns of U corresponding
  !          to indices of correct singular values should be correct.
  !          U is a two-dimensional REAL array, dimensioned U(NM,N).
  !
  !        V contains the matrix V (orthogonal) of the decomposition if
  !          MATV has been set to .TRUE.  Otherwise, V is not referenced.
  !          V may also coincide with A if U does not.  If an error
  !          exit is made, the columns of V corresponding to indices of
  !          correct singular values should be correct.  V is a two-
  !          dimensional REAL array, dimensioned V(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          K          if the K-th singular value has not been
  !                     determined after 30 iterations.
  !
  !        RV1 is a one-dimensional REAL array used for temporary storage,
  !          dimensioned RV1(N).
  !
  !     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***
  ! **See also:**  EISDOC
  !***
  ! **Routines called:**  PYTHAG

  !* REVISION HISTORY  (YYMMDD)
  !   811101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  !
  INTEGER i, j, k, l, M, N, ii, i1, kk, k1, ll, l1, mn, Nm, its, Ierr
  REAL A(Nm,*), W(*), U(Nm,*), V(Nm,*), Rv1(*)
  REAL c, f, g, h, s, x, y, z, scalee, s1
  LOGICAL Matu, Matv
  !
  !* FIRST EXECUTABLE STATEMENT  SVD
  Ierr = 0
  !
  DO i = 1, M
    !
    DO j = 1, N
      U(i,j) = A(i,j)
    END DO
  END DO
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
        scalee = scalee + ABS(U(k,i))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        !
        DO k = i, M
          U(k,i) = U(k,i)/scalee
          s = s + U(k,i)**2
        END DO
        !
        f = U(i,i)
        g = -SIGN(SQRT(s),f)
        h = f*g - s
        U(i,i) = f - g
        IF ( i/=N ) THEN
          !
          DO j = l, N
            s = 0.0E0
            !
            DO k = i, M
              s = s + U(k,i)*U(k,j)
            END DO
            !
            f = s/h
            !
            DO k = i, M
              U(k,j) = U(k,j) + f*U(k,i)
            END DO
          END DO
        END IF
        !
        DO k = i, M
          U(k,i) = scalee*U(k,i)
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
        scalee = scalee + ABS(U(i,k))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        !
        DO k = l, N
          U(i,k) = U(i,k)/scalee
          s = s + U(i,k)**2
        END DO
        !
        f = U(i,l)
        g = -SIGN(SQRT(s),f)
        h = f*g - s
        U(i,l) = f - g
        !
        DO k = l, N
          Rv1(k) = U(i,k)/h
        END DO
        !
        IF ( i/=M ) THEN
          !
          DO j = l, M
            s = 0.0E0
            !
            DO k = l, N
              s = s + U(j,k)*U(i,k)
            END DO
            !
            DO k = l, N
              U(j,k) = U(j,k) + s*Rv1(k)
            END DO
          END DO
        END IF
        !
        DO k = l, N
          U(i,k) = scalee*U(i,k)
        END DO
      END IF
    END IF
    !
    s1 = MAX(s1,ABS(W(i))+ABS(Rv1(i)))
  END DO
  !     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
  IF ( Matv ) THEN
    !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
    DO ii = 1, N
      i = N + 1 - ii
      IF ( i/=N ) THEN
        IF ( g/=0.0E0 ) THEN
          !
          DO j = l, N
            !     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            V(j,i) = (U(i,j)/U(i,l))/g
          END DO
          !
          DO j = l, N
            s = 0.0E0
            !
            DO k = l, N
              s = s + U(i,k)*V(k,j)
            END DO
            !
            DO k = l, N
              V(k,j) = V(k,j) + s*V(k,i)
            END DO
          END DO
        END IF
        !
        DO j = l, N
          V(i,j) = 0.0E0
          V(j,i) = 0.0E0
        END DO
      END IF
      !
      V(i,i) = 1.0E0
      g = Rv1(i)
      l = i
    END DO
  END IF
  !     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  IF ( Matu ) THEN
    !     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
    mn = N
    IF ( M<N ) mn = M
    !
    DO ii = 1, mn
      i = mn + 1 - ii
      l = i + 1
      g = W(i)
      IF ( i/=N ) THEN
        !
        DO j = l, N
          U(i,j) = 0.0E0
        END DO
      END IF
      !
      IF ( g==0.0E0 ) THEN
        !
        DO j = i, M
          U(j,i) = 0.0E0
        END DO
      ELSE
        IF ( i/=mn ) THEN
          !
          DO j = l, N
            s = 0.0E0
            !
            DO k = l, M
              s = s + U(k,i)*U(k,j)
            END DO
            !     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            f = (s/U(i,i))/g
            !
            DO k = i, M
              U(k,j) = U(k,j) + f*U(k,i)
            END DO
          END DO
        END IF
        !
        DO j = i, M
          U(j,i) = U(j,i)/g
          !
        END DO
      END IF
      !
      U(i,i) = U(i,i) + 1.0E0
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
      IF ( Matu ) THEN
        !
        DO j = 1, M
          y = U(j,l1)
          z = U(j,i)
          U(j,l1) = y*c + z*s
          U(j,i) = -y*s + z*c
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
        IF ( Matv ) THEN
          !
          DO j = 1, N
            x = V(j,i1)
            z = V(j,i)
            V(j,i1) = x*c + z*s
            V(j,i) = -x*s + z*c
          END DO
        END IF
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
        IF ( Matu ) THEN
          !
          DO j = 1, M
            y = U(j,i1)
            z = U(j,i)
            U(j,i1) = y*c + z*s
            U(j,i) = -y*s + z*c
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
      IF ( Matv ) THEN
        !
        DO j = 1, N
          V(j,k) = -V(j,k)
        END DO
      END IF
    END IF
    !
  END DO
  !
  RETURN
  !     .......... SET ERROR -- NO CONVERGENCE TO A
  !                SINGULAR VALUE AFTER 30 ITERATIONS ..........
  200  Ierr = k
  RETURN
END SUBROUTINE SVD
