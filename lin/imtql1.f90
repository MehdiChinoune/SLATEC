!*==IMTQL1.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK IMTQL1
      SUBROUTINE IMTQL1(N,D,E,Ierr)
      IMPLICIT NONE
!*--IMTQL15
!***BEGIN PROLOGUE  IMTQL1
!***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
!            using the implicit QL method.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (IMTQL1-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure IMTQL1,
!     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson,
!     as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
!
!     This subroutine finds the eigenvalues of a SYMMETRIC
!     TRIDIAGONAL matrix by the implicit QL method.
!
!     On INPUT
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!      On OUTPUT
!
!        D contains the eigenvalues in ascending order.  If an error
!          exit is made, the eigenvalues are correct and ordered for
!          indices 1, 2, ..., IERR-1, but may not be the smallest
!          eigenvalues.
!
!        E has been destroyed.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues should be correct for indices
!                     1, 2, ..., IERR-1.  These eigenvalues are
!                     ordered, but are not necessarily the smallest.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IMTQL1
!
      INTEGER i , j , l , m , N , ii , mml , Ierr
      REAL D(*) , E(*)
      REAL b , c , f , g , p , r , s , s1 , s2
      REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  IMTQL1
      Ierr = 0
      IF ( N/=1 ) THEN
!
        DO i = 2 , N
          E(i-1) = E(i)
        ENDDO
!
        E(N) = 0.0E0
!
        DO l = 1 , N
          j = 0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
 20       DO m = l , N
            IF ( m==N ) EXIT
            s1 = ABS(D(m)) + ABS(D(m+1))
            s2 = s1 + ABS(E(m))
            IF ( s2==s1 ) EXIT
          ENDDO
!
          p = D(l)
          IF ( m==l ) THEN
!     .......... ORDER EIGENVALUES ..........
            IF ( l/=1 ) THEN
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
              DO ii = 2 , l
                i = l + 2 - ii
                IF ( p>=D(i-1) ) GOTO 30
                D(i) = D(i-1)
              ENDDO
            ENDIF
!
            i = 1
 30         D(i) = p
          ELSE
            IF ( j==30 ) GOTO 100
            j = j + 1
!     .......... FORM SHIFT ..........
            g = (D(l+1)-p)/(2.0E0*E(l))
            r = PYTHAG(g,1.0E0)
            g = D(m) - p + E(l)/(g+SIGN(r,g))
            s = 1.0E0
            c = 1.0E0
            p = 0.0E0
            mml = m - l
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
            DO ii = 1 , mml
              i = m - ii
              f = s*E(i)
              b = c*E(i)
              IF ( ABS(f)<ABS(g) ) THEN
                s = f/g
                r = SQRT(s*s+1.0E0)
                E(i+1) = g*r
                c = 1.0E0/r
                s = s*c
              ELSE
                c = g/f
                r = SQRT(c*c+1.0E0)
                E(i+1) = f*r
                s = 1.0E0/r
                c = c*s
              ENDIF
              g = D(i+1) - p
              r = (D(i)-g)*s + 2.0E0*c*b
              p = s*r
              D(i+1) = g + p
              g = c*r - b
            ENDDO
!
            D(l) = D(l) - p
            E(l) = g
            E(m) = 0.0E0
            GOTO 20
          ENDIF
!
        ENDDO
      ENDIF
      GOTO 99999
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 100  Ierr = l
99999 END SUBROUTINE IMTQL1
