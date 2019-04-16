!** TQL2
SUBROUTINE TQL2(Nm,N,D,E,Z,Ierr)
  !>
  !***
  !  Compute the eigenvalues and eigenvectors of symmetric
  !            tridiagonal matrix.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A5, D4C2A
  !***
  ! **Type:**      SINGLE PRECISION (TQL2-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure TQL2,
  !     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
  !     Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
  !
  !     This subroutine finds the eigenvalues and eigenvectors
  !     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
  !     The eigenvectors of a FULL SYMMETRIC matrix can also
  !     be found if  TRED2  has been used to reduce this
  !     full matrix to tridiagonal form.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, Z, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is
  !          arbitrary.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        Z contains the transformation matrix produced in the
  !          reduction by  TRED2, if performed.  If the eigenvectors
  !          of the tridiagonal matrix are desired, Z must contain
  !          the identity matrix.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,N).
  !
  !      On Output
  !
  !        D contains the eigenvalues in ascending order.  If an
  !          error exit is made, the eigenvalues are correct but
  !          unordered for indices 1, 2, ..., IERR-1.
  !
  !        E has been destroyed.
  !
  !        Z contains orthonormal eigenvectors of the symmetric
  !          tridiagonal (or full) matrix.  If an error exit is made,
  !          Z contains the eigenvectors associated with the stored
  !          eigenvalues.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, m, N, ii, l1, l2, Nm, mml, Ierr
  REAL D(*), E(*), Z(Nm,*)
  REAL b, c, c2, c3, dl1, el1, f, g, h, p, r, s, s2
  !
  !* FIRST EXECUTABLE STATEMENT  TQL2
  Ierr = 0
  IF ( N/=1 ) THEN
    !
    DO i = 2, N
      E(i-1) = E(i)
    END DO
    !
    f = 0.0E0
    b = 0.0E0
    E(N) = 0.0E0
    !
    DO l = 1, N
      j = 0
      h = ABS(D(l)) + ABS(E(l))
      IF ( b<h ) b = h
      !     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
      DO m = l, N
        IF ( b+ABS(E(m))==b ) EXIT
        !     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
        !                THROUGH THE BOTTOM OF THE LOOP ..........
      END DO
      !
      IF ( m/=l ) THEN
        DO WHILE ( j/=30 )
          j = j + 1
          !     .......... FORM SHIFT ..........
          l1 = l + 1
          l2 = l1 + 1
          g = D(l)
          p = (D(l1)-g)/(2.0E0*E(l))
          r = PYTHAG(p,1.0E0)
          D(l) = E(l)/(p+SIGN(r,p))
          D(l1) = E(l)*(p+SIGN(r,p))
          dl1 = D(l1)
          h = g - D(l)
          IF ( l2<=N ) THEN
            !
            DO i = l2, N
              D(i) = D(i) - h
            END DO
          END IF
          !
          f = f + h
          !     .......... QL TRANSFORMATION ..........
          p = D(m)
          c = 1.0E0
          c2 = c
          el1 = E(l1)
          s = 0.0E0
          mml = m - l
          !     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
          DO ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c*E(i)
            h = c*p
            IF ( ABS(p)<ABS(E(i)) ) THEN
              c = p/E(i)
              r = SQRT(c*c+1.0E0)
              E(i+1) = s*E(i)*r
              s = 1.0E0/r
              c = c*s
            ELSE
              c = E(i)/p
              r = SQRT(c*c+1.0E0)
              E(i+1) = s*p*r
              s = c/r
              c = 1.0E0/r
            END IF
            p = c*D(i) - s*g
            D(i+1) = h + s*(c*g+s*D(i))
            !     .......... FORM VECTOR ..........
            DO k = 1, N
              h = Z(k,i+1)
              Z(k,i+1) = s*Z(k,i) + c*h
              Z(k,i) = c*Z(k,i) - s*h
            END DO
            !
          END DO
          !
          p = -s*s2*c3*el1*E(l)/dl1
          E(l) = s*p
          D(l) = c*p
          IF ( b+ABS(E(l))<=b ) GOTO 20
        END DO
        GOTO 100
      END IF
      20  D(l) = D(l) + f
    END DO
    !     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
    DO ii = 2, N
      i = ii - 1
      k = i
      p = D(i)
      !
      DO j = ii, N
        IF ( D(j)<p ) THEN
          k = j
          p = D(j)
        END IF
      END DO
      !
      IF ( k/=i ) THEN
        D(k) = D(i)
        D(i) = p
        !
        DO j = 1, N
          p = Z(j,i)
          Z(j,i) = Z(j,k)
          Z(j,k) = p
        END DO
      END IF
      !
      !
    END DO
  END IF
  RETURN
  !     .......... SET ERROR -- NO CONVERGENCE TO AN
  !                EIGENVALUE AFTER 30 ITERATIONS ..........
  100  Ierr = l
  RETURN
END SUBROUTINE TQL2
