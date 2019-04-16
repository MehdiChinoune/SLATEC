!** HTRIBK
SUBROUTINE HTRIBK(Nm,N,Ar,Ai,Tau,M,Zr,Zi)
  !>
  !***
  !  Form the eigenvectors of a complex Hermitian matrix from
  !            the eigenvectors of a real symmetric tridiagonal matrix
  !            output from HTRIDI.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      SINGLE PRECISION (HTRIBK-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of a complex analogue of
  !     the ALGOL procedure TRBAK1, NUM. MATH. 11, 181-195(1968)
  !     by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine forms the eigenvectors of a COMPLEX HERMITIAN
  !     matrix by back transforming those of the corresponding
  !     real symmetric tridiagonal matrix determined by  HTRIDI.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR, AI, ZR, and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        AR and AI contain some information about the unitary
  !          transformations used in the reduction by  HTRIDI  in the
  !          strict lower triangle of AR and the full lower triangle of
  !          AI.  The remaining upper parts of the matrices are arbitrary.
  !          AR and AI are two-dimensional REAL arrays, dimensioned
  !          AR(NM,N) and AI(NM,N).
  !
  !        TAU contains further information about the transformations.
  !          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
  !
  !        M is the number of eigenvectors to be back transformed.
  !          M is an INTEGER variable.
  !
  !       ZR contains the eigenvectors to be back transformed in its first
  !          M columns.  The contents of ZI are immaterial.  ZR and ZI are
  !          two-dimensional REAL arrays, dimensioned ZR(NM,M) and
  !          ZI(NM,M).
  !
  !     On OUTPUT
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the transformed eigenvectors in their first M columns.
  !
  !     Note that the last component of each returned vector
  !     is real and that vector Euclidean norms are preserved.
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
  INTEGER i, j, k, l, M, N, Nm
  REAL Ar(Nm,*), Ai(Nm,*), Tau(2,*), Zr(Nm,*), Zi(Nm,*)
  REAL h, s, si
  !
  !* FIRST EXECUTABLE STATEMENT  HTRIBK
  IF ( M/=0 ) THEN
    !     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
    !                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
    !                TRIDIAGONAL MATRIX. ..........
    DO k = 1, N
      !
      DO j = 1, M
        Zi(k,j) = -Zr(k,j)*Tau(2,k)
        Zr(k,j) = Zr(k,j)*Tau(1,k)
      END DO
    END DO
    !
    IF ( N/=1 ) THEN
      !     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO i = 2, N
        l = i - 1
        h = Ai(i,i)
        IF ( h/=0.0E0 ) THEN
          !
          DO j = 1, M
            s = 0.0E0
            si = 0.0E0
            !
            DO k = 1, l
              s = s + Ar(i,k)*Zr(k,j) - Ai(i,k)*Zi(k,j)
              si = si + Ar(i,k)*Zi(k,j) + Ai(i,k)*Zr(k,j)
            END DO
            !     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            s = (s/h)/h
            si = (si/h)/h
            !
            DO k = 1, l
              Zr(k,j) = Zr(k,j) - s*Ar(i,k) - si*Ai(i,k)
              Zi(k,j) = Zi(k,j) - si*Ar(i,k) + s*Ai(i,k)
            END DO
            !
          END DO
        END IF
        !
      END DO
    END IF
  END IF
  !
END SUBROUTINE HTRIBK
