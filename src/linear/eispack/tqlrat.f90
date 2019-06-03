!** TQLRAT
SUBROUTINE TQLRAT(N,D,E2,Ierr)
  !>
  !  Compute the eigenvalues of symmetric tridiagonal matrix
  !            using a rational variant of the QL method.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A5, D4C2A
  !***
  ! **Type:**      SINGLE PRECISION (TQLRAT-S)
  !***
  ! **Keywords:**  EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX, EISPACK,
  !             QL METHOD
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure TQLRAT.
  !
  !     This subroutine finds the eigenvalues of a SYMMETRIC
  !     TRIDIAGONAL matrix by the rational QL method.
  !
  !     On Input
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E2 contains the squares of the subdiagonal elements of the
  !          symmetric tridiagonal matrix in its last N-1 positions.
  !          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
  !          dimensioned E2(N).
  !
  !      On Output
  !
  !        D contains the eigenvalues in ascending order.  If an
  !          error exit is made, the eigenvalues are correct and
  !          ordered for indices 1, 2, ..., IERR-1, but may not be
  !          the smallest eigenvalues.
  !
  !        E2 has been destroyed.
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
  !               C. H. Reinsch, Eigenvalues of a REAL(SP), symmetric, tri-
  !                 diagonal matrix, Algorithm 464, Communications of the
  !                 ACM 16, 11 (November 1973), pp. 689.
  !***
  ! **Routines called:**  PYTHAG, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH
  INTEGER i, j, l, m, N, ii, l1, mml, Ierr
  REAL(SP) D(*), E2(*)
  REAL(SP) b, c, f, g, h, p, r, s
  !
  REAL(SP), PARAMETER :: machep = R1MACH(4)
  !* FIRST EXECUTABLE STATEMENT  TQLRAT
  !
  Ierr = 0
  IF ( N/=1 ) THEN
    !
    DO i = 2, N
      E2(i-1) = E2(i)
    END DO
    !
    f = 0.0E0
    b = 0.0E0
    E2(N) = 0.0E0
    !
    DO l = 1, N
      j = 0
      h = machep*(ABS(D(l))+SQRT(E2(l)))
      IF ( b<=h ) THEN
        b = h
        c = b*b
      END IF
      !     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
      DO m = l, N
        IF ( E2(m)<=c ) EXIT
        !     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
        !                THROUGH THE BOTTOM OF THE LOOP ..........
      END DO
      !
      IF ( m/=l ) THEN
        DO WHILE ( j/=30 )
          j = j + 1
          !     .......... FORM SHIFT ..........
          l1 = l + 1
          s = SQRT(E2(l))
          g = D(l)
          p = (D(l1)-g)/(2.0E0*s)
          r = PYTHAG(p,1.0E0)
          D(l) = s/(p+SIGN(r,p))
          h = g - D(l)
          !
          DO i = l1, N
            D(i) = D(i) - h
          END DO
          !
          f = f + h
          !     .......... RATIONAL QL TRANSFORMATION ..........
          g = D(m)
          IF ( g==0.0E0 ) g = b
          h = g
          s = 0.0E0
          mml = m - l
          !     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
          DO ii = 1, mml
            i = m - ii
            p = g*h
            r = p + E2(i)
            E2(i+1) = s*r
            s = E2(i)/r
            D(i+1) = h + s*(h+D(i))
            g = D(i) - E2(i)/g
            IF ( g==0.0E0 ) g = b
            h = g*p/r
          END DO
          !
          E2(l) = s*g
          D(l) = h
          !     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
          IF ( h==0.0E0 ) GOTO 20
          IF ( ABS(E2(l))<=ABS(c/h) ) GOTO 20
          E2(l) = h*E2(l)
          IF ( E2(l)==0.0E0 ) GOTO 20
        END DO
        GOTO 100
      END IF
      20  p = D(l) + f
      !     .......... ORDER EIGENVALUES ..........
      IF ( l/=1 ) THEN
        !     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
        DO ii = 2, l
          i = l + 2 - ii
          IF ( p>=D(i-1) ) GOTO 40
          D(i) = D(i-1)
        END DO
      END IF
      !
      i = 1
      40  D(i) = p
      !
    END DO
  END IF
  RETURN
  !     .......... SET ERROR -- NO CONVERGENCE TO AN
  !                EIGENVALUE AFTER 30 ITERATIONS ..........
  100  Ierr = l
  RETURN
END SUBROUTINE TQLRAT
