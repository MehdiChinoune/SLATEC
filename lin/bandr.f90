!DECK BANDR
SUBROUTINE BANDR(Nm,N,Mb,A,D,E,E2,Matz,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BANDR
  !***PURPOSE  Reduce a real symmetric band matrix to symmetric
  !            tridiagonal matrix and, optionally, accumulate
  !            orthogonal similarity transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B1
  !***TYPE      SINGLE PRECISION (BANDR-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure BANDRD,
  !     NUM. MATH. 12, 231-241(1968) by Schwarz.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 273-283(1971).
  !
  !     This subroutine reduces a REAL SYMMETRIC BAND matrix
  !     to a symmetric tridiagonal matrix using and optionally
  !     accumulating orthogonal similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        MB is the (half) band width of the matrix, defined as the
  !          number of adjacent diagonals, including the principal
  !          diagonal, required to specify the non-zero portion of the
  !          lower triangle of the matrix.  MB is less than or equal
  !          to N.  MB is an INTEGER variable.
  !
  !        A contains the lower triangle of the real symmetric band
  !          matrix.  Its lowest subdiagonal is stored in the last
  !          N+1-MB  positions of the first column, its next subdiagonal
  !          in the last  N+2-MB  positions of the second column, further
  !          subdiagonals similarly, and finally its principal diagonal
  !          in the  N  positions of the last column.  Contents of storage
  !          locations not part of the matrix are arbitrary.  A is a
  !          two-dimensional REAL array, dimensioned A(NM,MB).
  !
  !        MATZ should be set to .TRUE. if the transformation matrix is
  !          to be accumulated, and to .FALSE. otherwise.  MATZ is a
  !          LOGICAL variable.
  !
  !     On OUTPUT
  !
  !        A has been destroyed, except for its last two columns which
  !          contain a copy of the tridiagonal matrix.
  !
  !        D contains the diagonal elements of the tridiagonal matrix.
  !          D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the tridiagonal
  !          matrix in its last N-1 positions.  E(1) is set to zero.
  !          E is a one-dimensional REAL array, dimensioned E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2 may coincide with E if the squares are not needed.
  !          E2 is a one-dimensional REAL array, dimensioned E2(N).
  !
  !        Z contains the orthogonal transformation matrix produced in
  !          the reduction if MATZ has been set to .TRUE.  Otherwise, Z
  !          is not referenced.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,N).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  BANDR
  !
  INTEGER j, k, l, N, r, i1, i2, j1, j2, kr, Mb, mr, m1, Nm, &
    n2, r1, ugl, maxl, maxr
  REAL A(Nm,*), D(*), E(*), E2(*), Z(Nm,*)
  REAL g, u, b1, b2, c2, f1, f2, s2, dmin, dminrt
  LOGICAL Matz
  !
  !***FIRST EXECUTABLE STATEMENT  BANDR
  dmin = 2.0E0**(-64)
  dminrt = 2.0E0**(-32)
  !     .......... INITIALIZE DIAGONAL SCALING MATRIX ..........
  DO j = 1, N
    D(j) = 1.0E0
  ENDDO
  !
  IF ( Matz ) THEN
    !
    DO j = 1, N
      !
      DO k = 1, N
        Z(j,k) = 0.0E0
      ENDDO
      !
      Z(j,j) = 1.0E0
    ENDDO
  ENDIF
  !
  m1 = Mb - 1
  IF ( m1<1 ) THEN
    !
    DO j = 1, N
      D(j) = A(j,Mb)
      E(j) = 0.0E0
      E2(j) = 0.0E0
    ENDDO
    GOTO 99999
  ELSEIF ( m1/=1 ) THEN
    n2 = N - 2
    !
    DO k = 1, n2
      maxr = MIN(m1,N-k)
      !     .......... FOR R=MAXR STEP -1 UNTIL 2 DO -- ..........
      DO r1 = 2, maxr
        r = maxr + 2 - r1
        kr = k + r
        mr = Mb - r
        g = A(kr,mr)
        A(kr-1,1) = A(kr-1,mr+1)
        ugl = k
        !
        DO j = kr, N, m1
          j1 = j - 1
          j2 = j1 - 1
          IF ( g==0.0E0 ) EXIT
          b1 = A(j1,1)/g
          b2 = b1*D(j1)/D(j)
          s2 = 1.0E0/(1.0E0+b1*b2)
          IF ( s2>=0.5E0 ) THEN
            !
            u = D(j1)
            D(j1) = s2*D(j)
            D(j) = s2*u
            f1 = 2.0E0*A(j,m1)
            f2 = b1*A(j,Mb)
            u = b1*(f2-f1) + A(j1,Mb)
            A(j,m1) = b2*(b1*A(j,m1)-A(j1,Mb)) + f2 - A(j,m1)
            A(j1,Mb) = b2*(b2*A(j1,Mb)+f1) + A(j,Mb)
            A(j,Mb) = u
            !
            DO l = ugl, j2
              i2 = Mb - j + l
              u = b2*A(j1,i2+1) + A(j,i2)
              A(j,i2) = -A(j1,i2+1) + b1*A(j,i2)
              A(j1,i2+1) = u
            ENDDO
            !
            ugl = j
            A(j1,1) = b2*A(j1,1) + g
            IF ( j/=N ) THEN
              maxl = MIN(m1,N-j1)
              !
              DO l = 2, maxl
                i1 = j1 + l
                i2 = Mb - l
                u = b2*A(i1,i2) + A(i1,i2+1)
                A(i1,i2+1) = -A(i1,i2) + b1*A(i1,i2+1)
                A(i1,i2) = u
              ENDDO
              !
              i1 = j + m1
              IF ( i1<=N ) THEN
                g = A(i1,1)
                A(i1,1) = b1*A(i1,1)
              ENDIF
            ENDIF
            IF ( Matz ) THEN
              !
              DO l = 1, N
                u = b2*Z(l,j1) + Z(l,j)
                Z(l,j) = -Z(l,j1) + b1*Z(l,j)
                Z(l,j1) = u
              ENDDO
            ENDIF
          ELSE
            b1 = g/A(j1,1)
            b2 = b1*D(j)/D(j1)
            c2 = 1.0E0 - s2
            D(j1) = c2*D(j1)
            D(j) = c2*D(j)
            f1 = 2.0E0*A(j,m1)
            f2 = b1*A(j1,Mb)
            A(j,m1) = -b2*(b1*A(j,m1)-A(j,Mb)) - f2 + A(j,m1)
            A(j1,Mb) = b2*(b2*A(j,Mb)+f1) + A(j1,Mb)
            A(j,Mb) = b1*(f2-f1) + A(j,Mb)
            !
            DO l = ugl, j2
              i2 = Mb - j + l
              u = A(j1,i2+1) + b2*A(j,i2)
              A(j,i2) = -b1*A(j1,i2+1) + A(j,i2)
              A(j1,i2+1) = u
            ENDDO
            !
            ugl = j
            A(j1,1) = A(j1,1) + b2*g
            IF ( j/=N ) THEN
              maxl = MIN(m1,N-j1)
              !
              DO l = 2, maxl
                i1 = j1 + l
                i2 = Mb - l
                u = A(i1,i2) + b2*A(i1,i2+1)
                A(i1,i2+1) = -b1*A(i1,i2) + A(i1,i2+1)
                A(i1,i2) = u
              ENDDO
              !
              i1 = j + m1
              IF ( i1<=N ) g = b2*A(i1,1)
            ENDIF
            IF ( Matz ) THEN
              !
              DO l = 1, N
                u = Z(l,j1) + b2*Z(l,j)
                Z(l,j) = -b1*Z(l,j1) + Z(l,j)
                Z(l,j1) = u
                !
              ENDDO
            ENDIF
          ENDIF
          !
        ENDDO
        !
      ENDDO
      !
      IF ( MOD(k,64)==0 ) THEN
        !     .......... RESCALE TO AVOID UNDERFLOW OR OVERFLOW ..........
        DO j = k, N
          IF ( D(j)<dmin ) THEN
            maxl = MAX(1,Mb+1-j)
            !
            DO l = maxl, m1
              A(j,l) = dminrt*A(j,l)
            ENDDO
            !
            IF ( j/=N ) THEN
              maxl = MIN(m1,N-j)
              !
              DO l = 1, maxl
                i1 = j + l
                i2 = Mb - l
                A(i1,i2) = dminrt*A(i1,i2)
              ENDDO
            ENDIF
            !
            IF ( Matz ) THEN
              !
              DO l = 1, N
                Z(l,j) = dminrt*Z(l,j)
              ENDDO
            ENDIF
            !
            A(j,Mb) = dmin*A(j,Mb)
            D(j) = D(j)/dmin
          ENDIF
        ENDDO
      ENDIF
      !
    ENDDO
  ENDIF
  !     .......... FORM SQUARE ROOT OF SCALING MATRIX ..........
  DO j = 2, N
    E(j) = SQRT(D(j))
  ENDDO
  !
  IF ( Matz ) THEN
    !
    DO j = 1, N
      !
      DO k = 2, N
        Z(j,k) = E(k)*Z(j,k)
      ENDDO
      !
    ENDDO
  ENDIF
  !
  u = 1.0E0
  !
  DO j = 2, N
    A(j,m1) = u*E(j)*A(j,m1)
    u = E(j)
    E2(j) = A(j,m1)**2
    A(j,Mb) = D(j)*A(j,Mb)
    D(j) = A(j,Mb)
    E(j) = A(j,m1)
  ENDDO
  !
  D(1) = A(1,Mb)
  E(1) = 0.0E0
  E2(1) = 0.0E0
  !
  99999 CONTINUE
  END SUBROUTINE BANDR
