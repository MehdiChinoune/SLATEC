!*==BQR.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK BQR
      SUBROUTINE BQR(Nm,N,Mb,A,T,R,Ierr,Nv,Rv)
      IMPLICIT NONE
!*--BQR5
!***BEGIN PROLOGUE  BQR
!***PURPOSE  Compute some of the eigenvalues of a real symmetric
!            matrix using the QR method with shifts of origin.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A6
!***TYPE      SINGLE PRECISION (BQR-S)
!***KEYWORDS  EIGENVALUES, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure BQR,
!     NUM. MATH. 16, 85-92(1970) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 266-272(1971).
!
!     This subroutine finds the eigenvalue of smallest (usually)
!     magnitude of a REAL SYMMETRIC BAND matrix using the
!     QR algorithm with shifts of origin.  Consecutive calls
!     can be made to find further eigenvalues.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        MB is the (half) band width of the matrix, defined as the
!          number of adjacent diagonals, including the principal
!          diagonal, required to specify the non-zero portion of the
!          lower triangle of the matrix.  MB is an INTEGER variable.
!          MB must be less than or equal to N on first call.
!
!        A contains the lower triangle of the symmetric band input
!          matrix stored as an N by MB array.  Its lowest subdiagonal
!          is stored in the last N+1-MB positions of the first column,
!          its next subdiagonal in the last N+2-MB positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the N positions of the last column.
!          Contents of storages not part of the matrix are arbitrary.
!          On a subsequent call, its output contents from the previous
!          call should be passed.  A is a two-dimensional REAL array,
!          dimensioned A(NM,MB).
!
!        T specifies the shift (of eigenvalues) applied to the diagonal
!          of A in forming the input matrix. What is actually determined
!          is the eigenvalue of A+TI (I is the identity matrix) nearest
!          to T.  On a subsequent call, the output value of T from the
!          previous call should be passed if the next nearest eigenvalue
!          is sought.  T is a REAL variable.
!
!        R should be specified as zero on the first call, and as its
!          output value from the previous call on a subsequent call.
!          It is used to determine when the last row and column of
!          the transformed band matrix can be regarded as negligible.
!          R is a REAL variable.
!
!        NV must be set to the dimension of the array parameter RV
!          as declared in the calling program dimension statement.
!          NV is an INTEGER variable.
!
!     On OUTPUT
!
!        A contains the transformed band matrix.  The matrix A+TI
!          derived from the output parameters is similar to the
!          input A+TI to within rounding errors.  Its last row and
!          column are null (if IERR is zero).
!
!        T contains the computed eigenvalue of A+TI (if IERR is zero),
!          where I is the identity matrix.
!
!        R contains the maximum of its input value and the norm of the
!          last column of the input matrix A.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30 iterations.
!
!        RV is a one-dimensional REAL array of dimension NV which is
!          at least (2*MB**2+4*MB-3), used for temporary storage.  The
!          first (3*MB-2) locations correspond to the ALGOL array B,
!          the next (2*MB-1) locations correspond to the ALGOL array H,
!          and the final (2*MB**2-MB) locations correspond to the MB
!          by (2*MB-1) ALGOL array U.
!
!     NOTE. For a subsequent call, N should be replaced by N-1, but
!     MB should not be altered even when it exceeds the current N.
!
!     Calls PYTHAG(A,B) for SQRT(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BQR
!
      INTEGER i , j , k , l , m , N , ii , ik , jk , jm , kj , kk , km , ll , 
     &        Mb , mk , mn , mz
      INTEGER m1 , m2 , m3 , m4 , ni , Nm , Nv , its , kj1 , m21 , m31 , Ierr , 
     &        imult
      REAL A(Nm,*) , Rv(*)
      REAL f , g , q , R , s , T , scale
      REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  BQR
      Ierr = 0
      m1 = MIN(Mb,N)
      m = m1 - 1
      m2 = m + m
      m21 = m2 + 1
      m3 = m21 + m
      m31 = m3 + 1
      m4 = m31 + m2
      mn = m + N
      mz = Mb - m1
      its = 0
      DO
!     .......... TEST FOR CONVERGENCE ..........
        g = A(N,Mb)
        IF ( m==0 ) EXIT
        f = 0.0E0
!
        DO k = 1 , m
          mk = k + mz
          f = f + ABS(A(N,mk))
        ENDDO
!
        IF ( its==0.AND.f>R ) R = f
        IF ( R+f<=R ) EXIT
        IF ( its==30 ) THEN
!     .......... SET ERROR -- NO CONVERGENCE TO
!                EIGENVALUE AFTER 30 ITERATIONS ..........
          Ierr = N
          GOTO 99999
        ELSE
          its = its + 1
!     .......... FORM SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
          IF ( f<=0.25E0*R.OR.its>=5 ) THEN
            f = A(N,Mb-1)
            IF ( f/=0.0E0 ) THEN
              q = (A(N-1,Mb)-g)/(2.0E0*f)
              s = PYTHAG(q,1.0E0)
              g = g - f/(q+SIGN(s,q))
            ENDIF
            T = T + g
!
            DO i = 1 , N
              A(i,Mb) = A(i,Mb) - g
            ENDDO
          ENDIF
!
          DO k = m31 , m4
            Rv(k) = 0.0E0
          ENDDO
!
          DO ii = 1 , mn
            i = ii - m
            ni = N - ii
            IF ( ni<0 ) GOTO 20
!     .......... FORM COLUMN OF SHIFTED MATRIX A-G*I ..........
            l = MAX(1,2-i)
!
            DO k = 1 , m3
              Rv(k) = 0.0E0
            ENDDO
!
            DO k = l , m1
              km = k + m
              mk = k + mz
              Rv(km) = A(ii,mk)
            ENDDO
!
            ll = MIN(m,ni)
            IF ( ll/=0 ) THEN
!
              DO k = 1 , ll
                km = k + m21
                ik = ii + k
                mk = Mb - k
                Rv(km) = A(ik,mk)
              ENDDO
            ENDIF
!     .......... PRE-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
            ll = m2
            imult = 0
!     .......... MULTIPLICATION PROCEDURE ..........
 10         kj = m4 - m1
!
            DO j = 1 , ll
              kj = kj + m1
              jm = j + m3
              IF ( Rv(jm)/=0.0E0 ) THEN
                f = 0.0E0
!
                DO k = 1 , m1
                  kj = kj + 1
                  jk = j + k - 1
                  f = f + Rv(kj)*Rv(jk)
                ENDDO
!
                f = f/Rv(jm)
                kj = kj - m1
!
                DO k = 1 , m1
                  kj = kj + 1
                  jk = j + k - 1
                  Rv(jk) = Rv(jk) - Rv(kj)*f
                ENDDO
!
                kj = kj - m1
              ENDIF
            ENDDO
!
            IF ( imult/=0 ) THEN
!     .......... STORE COLUMN OF NEW A MATRIX ..........
              DO k = l , m1
                mk = k + mz
                A(i,mk) = Rv(k)
              ENDDO
              GOTO 30
            ELSE
!     .......... HOUSEHOLDER REFLECTION ..........
              f = Rv(m21)
              s = 0.0E0
              Rv(m4) = 0.0E0
              scale = 0.0E0
!
              DO k = m21 , m3
                scale = scale + ABS(Rv(k))
              ENDDO
!
              IF ( scale/=0.0E0 ) THEN
!
                DO k = m21 , m3
                  s = s + (Rv(k)/scale)**2
                ENDDO
!
                s = scale*scale*s
                g = -SIGN(SQRT(s),f)
                Rv(m21) = g
                Rv(m4) = s - f*g
                kj = m4 + m2*m1 + 1
                Rv(kj) = f - g
!
                DO k = 2 , m1
                  kj = kj + 1
                  km = k + m2
                  Rv(kj) = Rv(km)
                ENDDO
              ENDIF
!     .......... SAVE COLUMN OF TRIANGULAR FACTOR R ..........
              DO k = l , m1
                km = k + m
                mk = k + mz
                A(ii,mk) = Rv(km)
              ENDDO
            ENDIF
!
 20         l = MAX(1,m1+1-i)
            IF ( i>0 ) THEN
!     .......... PERFORM ADDITIONAL STEPS ..........
              DO k = 1 , m21
                Rv(k) = 0.0E0
              ENDDO
!
              ll = MIN(m1,ni+m1)
!     .......... GET ROW OF TRIANGULAR FACTOR R ..........
              DO kk = 1 , ll
                k = kk - 1
                km = k + m1
                ik = i + k
                mk = Mb - k
                Rv(km) = A(ik,mk)
              ENDDO
!     .......... POST-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
              ll = m1
              imult = 1
              GOTO 10
            ENDIF
!     .......... UPDATE HOUSEHOLDER REFLECTIONS ..........
 30         IF ( l>1 ) l = l - 1
            kj1 = m4 + l*m1
!
            DO j = l , m2
              jm = j + m3
              Rv(jm) = Rv(jm+1)
!
              DO k = 1 , m1
                kj1 = kj1 + 1
                kj = kj1 - m1
                Rv(kj) = Rv(kj1)
              ENDDO
            ENDDO
!
!
          ENDDO
        ENDIF
      ENDDO
!     .......... CONVERGENCE ..........
      T = T + g
!
      DO i = 1 , N
        A(i,Mb) = A(i,Mb) - g
      ENDDO
!
      DO k = 1 , m1
        mk = k + mz
        A(N,mk) = 0.0E0
!
      ENDDO
99999 END SUBROUTINE BQR
