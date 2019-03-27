!** BANDV
SUBROUTINE BANDV(Nm,N,Mbw,A,E21,M,W,Z,Ierr,Nv,Rv,Rv6)
  IMPLICIT NONE
  !>
  !***
  !  Form the eigenvectors of a real symmetric band matrix
  !            associated with a set of ordered approximate eigenvalues
  !            by inverse iteration.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C3
  !***
  ! **Type:**      SINGLE PRECISION (BANDV-S)
  !***
  ! **Keywords:**  EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine finds those eigenvectors of a REAL SYMMETRIC
  !     BAND matrix corresponding to specified eigenvalues, using inverse
  !     iteration.  The subroutine may also be used to solve systems
  !     of linear equations with a symmetric or non-symmetric band
  !     coefficient matrix.
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
  !        MBW is the number of columns of the array A used to store the
  !          band matrix.  If the matrix is symmetric, MBW is its (half)
  !          band width, denoted MB and defined as the number of adjacent
  !          diagonals, including the principal diagonal, required to
  !          specify the non-zero portion of the lower triangle of the
  !          matrix.  If the subroutine is being used to solve systems
  !          of linear equations and the coefficient matrix is not
  !          symmetric, it must however have the same number of adjacent
  !          diagonals above the main diagonal as below, and in this
  !          case, MBW=2*MB-1.  MBW is an INTEGER variable.  MB must not
  !          be greater than N.
  !
  !        A contains the lower triangle of the symmetric band input
  !          matrix stored as an N by MB array.  Its lowest subdiagonal
  !          is stored in the last N+1-MB positions of the first column,
  !          its next subdiagonal in the last N+2-MB positions of the
  !          second column, further subdiagonals similarly, and finally
  !          its principal diagonal in the N positions of column MB.
  !          If the subroutine is being used to solve systems of linear
  !          equations and the coefficient matrix is not symmetric, A is
  !          N by 2*MB-1 instead with lower triangle as above and with
  !          its first superdiagonal stored in the first N-1 positions of
  !          column MB+1, its second superdiagonal in the first N-2
  !          positions of column MB+2, further superdiagonals similarly,
  !          and finally its highest superdiagonal in the first N+1-MB
  !          positions of the last column.  Contents of storage locations
  !          not part of the matrix are arbitrary.  A is a two-dimensional
  !          REAL array, dimensioned A(NM,MBW).
  !
  !        E21 specifies the ordering of the eigenvalues and contains
  !            0.0E0 if the eigenvalues are in ascending order, or
  !            2.0E0 if the eigenvalues are in descending order.
  !          If the subroutine is being used to solve systems of linear
  !          equations, E21 should be set to 1.0E0 if the coefficient
  !          matrix is symmetric and to -1.0E0 if not.  E21 is a REAL
  !          variable.
  !
  !        M is the number of specified eigenvalues or the number of
  !          systems of linear equations.  M is an INTEGER variable.
  !
  !        W contains the M eigenvalues in ascending or descending order.
  !          If the subroutine is being used to solve systems of linear
  !          equations (A-W(J)*I)*X(J)=B(J), where I is the identity
  !          matrix, W(J) should be set accordingly, for J=1,2,...,M.
  !          W is a one-dimensional REAL array, dimensioned W(M).
  !
  !        Z contains the constant matrix columns (B(J),J=1,2,...,M), if
  !          the subroutine is used to solve systems of linear equations.
  !          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
  !
  !        NV must be set to the dimension of the array parameter RV
  !          as declared in the calling program dimension statement.
  !          NV is an INTEGER variable.
  !
  !     On OUTPUT
  !
  !        A and W are unaltered.
  !
  !        Z contains the associated set of orthogonal eigenvectors.
  !          Any vector which fails to converge is set to zero.  If the
  !          subroutine is used to solve systems of linear equations,
  !          Z contains the solution matrix columns (X(J),J=1,2,...,M).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          -J         if the eigenvector corresponding to the J-th
  !                     eigenvalue fails to converge, or if the J-th
  !                     system of linear equations is nearly singular.
  !
  !        RV and RV6 are temporary storage arrays.  If the subroutine
  !          is being used to solve systems of linear equations, the
  !          determinant (up to sign) of A-W(M)*I is available, upon
  !          return, as the product of the first N elements of RV.
  !          RV and RV6 are one-dimensional REAL arrays.  Note that RV
  !          is dimensioned RV(NV), where NV must be at least N*(2*MB-1).
  !          RV6 is dimensioned RV6(N).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  !
  INTEGER i, j, k, M, N, r, ii, ij, jj, kj, mb, m1, Nm, Nv, &
    ij1, its, kj1, Mbw, m21
  INTEGER Ierr, maxj, maxk, group
  REAL A(Nm,*), W(*), Z(Nm,*), Rv(*), Rv6(*)
  REAL u, v, uk, xu, x0, x1, E21, eps2, eps3, eps4, norm, order, s
  !
  !* FIRST EXECUTABLE STATEMENT  BANDV
  Ierr = 0
  IF ( M/=0 ) THEN
    mb = Mbw
    IF ( E21<0.0E0 ) mb = (Mbw+1)/2
    m1 = mb - 1
    m21 = m1 + mb
    order = 1.0E0 - ABS(E21)
    !     .......... FIND VECTORS BY INVERSE ITERATION ..........
    DO r = 1, M
      its = 1
      x1 = W(r)
      IF ( r==1 ) THEN
        !     .......... COMPUTE NORM OF MATRIX ..........
        norm = 0.0E0
        !
        DO j = 1, mb
          jj = mb + 1 - j
          kj = jj + m1
          ij = 1
          s = 0.0E0
          !
          DO i = jj, N
            s = s + ABS(A(i,j))
            IF ( E21<0.0E0 ) THEN
              s = s + ABS(A(ij,kj))
              ij = ij + 1
            ENDIF
          ENDDO
          !
          norm = MAX(norm,s)
        ENDDO
        !
        IF ( E21<0.0E0 ) norm = 0.5E0*norm
        !     .......... EPS2 IS THE CRITERION FOR GROUPING,
        !                EPS3 REPLACES ZERO PIVOTS AND EQUAL
        !                ROOTS ARE MODIFIED BY EPS3,
        !                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
        IF ( norm==0.0E0 ) norm = 1.0E0
        eps2 = 1.0E-3*norm*ABS(order)
        eps3 = norm
        DO
          eps3 = 0.5E0*eps3
          IF ( norm+eps3<=norm ) THEN
            uk = SQRT(REAL(N))
            eps3 = uk*eps3
            eps4 = uk*eps3
            group = 0
            EXIT
          ENDIF
        ENDDO
        !     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
      ELSEIF ( ABS(x1-x0)>=eps2 ) THEN
        group = 0
      ELSE
        group = group + 1
        IF ( order*(x1-x0)<=0.0E0 ) x1 = x0 + order*eps3
      ENDIF
      !     .......... EXPAND MATRIX, SUBTRACT EIGENVALUE,
      !                AND INITIALIZE VECTOR ..........
      DO i = 1, N
        ij = i + MIN(0,i-m1)*N
        kj = ij + mb*N
        ij1 = kj + m1*N
        IF ( m1/=0 ) THEN
          !
          DO j = 1, m1
            IF ( ij>m1 ) THEN
              Rv(ij) = A(i,j)
            ELSEIF ( ij<=0 ) THEN
              Rv(ij1) = 0.0E0
              ij1 = ij1 + N
            ENDIF
            ij = ij + N
            ii = i + j
            IF ( ii<=N ) THEN
              jj = mb - j
              IF ( E21<0.0E0 ) THEN
                ii = i
                jj = mb + j
              ENDIF
              Rv(kj) = A(ii,jj)
              kj = kj + N
            ENDIF
          ENDDO
        ENDIF
        !
        Rv(ij) = A(i,mb) - x1
        Rv6(i) = eps4
        IF ( order==0.0E0 ) Rv6(i) = Z(i,r)
      ENDDO
      !
      IF ( m1/=0 ) THEN
        !     .......... ELIMINATION WITH INTERCHANGES ..........
        DO i = 1, N
          ii = i + 1
          maxk = MIN(i+m1-1,N)
          maxj = MIN(N-i,m21-2)*N
          !
          DO k = i, maxk
            kj1 = k
            j = kj1 + N
            jj = j + maxj
            !
            DO kj = j, jj, N
              Rv(kj1) = Rv(kj)
              kj1 = kj
            ENDDO
            !
            Rv(kj1) = 0.0E0
          ENDDO
          !
          IF ( i/=N ) THEN
            u = 0.0E0
            maxk = MIN(i+m1,N)
            maxj = MIN(N-ii,m21-2)*N
            !
            DO j = i, maxk
              IF ( ABS(Rv(j))>=ABS(u) ) THEN
                u = Rv(j)
                k = j
              ENDIF
            ENDDO
            !
            j = i + N
            jj = j + maxj
            IF ( k/=i ) THEN
              kj = k
              !
              DO ij = i, jj, N
                v = Rv(ij)
                Rv(ij) = Rv(kj)
                Rv(kj) = v
                kj = kj + N
              ENDDO
              !
              IF ( order==0.0E0 ) THEN
                v = Rv6(i)
                Rv6(i) = Rv6(k)
                Rv6(k) = v
              ENDIF
            ENDIF
            IF ( u/=0.0E0 ) THEN
              !
              DO k = ii, maxk
                v = Rv(k)/u
                kj = k
                !
                DO ij = j, jj, N
                  kj = kj + N
                  Rv(kj) = Rv(kj) - v*Rv(ij)
                ENDDO
                !
                IF ( order==0.0E0 ) Rv6(k) = Rv6(k) - v*Rv6(i)
              ENDDO
            ENDIF
          ENDIF
          !
        ENDDO
      ENDIF
      DO
        !     .......... BACK SUBSTITUTION
        !                FOR I=N STEP -1 UNTIL 1 DO -- ..........
        DO ii = 1, N
          i = N + 1 - ii
          maxj = MIN(ii,m21)
          IF ( maxj/=1 ) THEN
            ij1 = i
            j = ij1 + N
            jj = j + (maxj-2)*N
            !
            DO ij = j, jj, N
              ij1 = ij1 + 1
              Rv6(i) = Rv6(i) - Rv(ij)*Rv6(ij1)
            ENDDO
          ENDIF
          !
          v = Rv(i)
          IF ( ABS(v)<eps3 ) THEN
            !     .......... SET ERROR -- NEARLY SINGULAR LINEAR SYSTEM ..........
            IF ( order==0.0E0 ) Ierr = -r
            v = SIGN(eps3,v)
          ENDIF
          Rv6(i) = Rv6(i)/v
        ENDDO
        !
        xu = 1.0E0
        IF ( order==0.0E0 ) EXIT
        !     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
        !                MEMBERS OF GROUP ..........
        IF ( group/=0 ) THEN
          !
          DO jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0E0
            !
            DO i = 1, N
              xu = xu + Rv6(i)*Z(i,j)
            ENDDO
            !
            DO i = 1, N
              Rv6(i) = Rv6(i) - xu*Z(i,j)
            ENDDO
            !
          ENDDO
        ENDIF
        !
        norm = 0.0E0
        !
        DO i = 1, N
          norm = norm + ABS(Rv6(i))
        ENDDO
        !
        IF ( norm>=0.1E0 ) THEN
          !     .......... NORMALIZE SO THAT SUM OF SQUARES IS
          !                1 AND EXPAND TO FULL ORDER ..........
          u = 0.0E0
          !
          DO i = 1, N
            u = u + Rv6(i)**2
          ENDDO
          !
          xu = 1.0E0/SQRT(u)
          EXIT
          !     .......... IN-LINE PROCEDURE FOR CHOOSING
          !                A NEW STARTING VECTOR ..........
        ELSEIF ( its>=N ) THEN
          !     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
          Ierr = -r
          xu = 0.0E0
          EXIT
        ELSE
          its = its + 1
          xu = eps4/(uk+1.0E0)
          Rv6(1) = eps4
          !
          DO i = 2, N
            Rv6(i) = xu
          ENDDO
          !
          Rv6(its) = Rv6(its) - eps4*uk
        ENDIF
      ENDDO
      !
      DO i = 1, N
        Z(i,r) = Rv6(i)*xu
      ENDDO
      !
      x0 = x1
    ENDDO
  ENDIF
  !
END SUBROUTINE BANDV
