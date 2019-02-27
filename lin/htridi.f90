!DECK HTRIDI
SUBROUTINE HTRIDI(Nm,N,Ar,Ai,D,E,E2,Tau)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  HTRIDI
  !***PURPOSE  Reduce a complex Hermitian matrix to a real symmetric
  !            tridiagonal matrix using unitary similarity
  !            transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B1
  !***TYPE      SINGLE PRECISION (HTRIDI-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of a complex analogue of
  !     the ALGOL procedure TRED1, NUM. MATH. 11, 181-195(1968)
  !     by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine reduces a COMPLEX HERMITIAN matrix
  !     to a real symmetric tridiagonal matrix using
  !     unitary similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR and AI, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A=(AR,AI).  N is an INTEGER
  !          variable. N must be less than or equal to NM.
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the complex Hermitian input matrix.  Only the lower
  !          triangle of the matrix need be supplied.  AR and AI are two-
  !          dimensional REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
  !
  !     On OUTPUT
  !
  !        AR and AI contain some information about the unitary trans-
  !          formations used in the reduction in the strict lower triangle
  !          of AR and the full lower triangle of AI.  The rest of the
  !          matrices are unaltered.
  !
  !        D contains the diagonal elements of the real symmetric
  !          tridiagonal matrix.  D is a one-dimensional REAL array,
  !          dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the real tridiagonal
  !          matrix in its last N-1 positions.  E(1) is set to zero.
  !          E is a one-dimensional REAL array, dimensioned E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2(1) is set to zero.  E2 may coincide with E if the squares
  !          are not needed.  E2 is a one-dimensional REAL array,
  !          dimensioned E2(N).
  !
  !        TAU contains further information about the transformations.
  !          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
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
  !***END PROLOGUE  HTRIDI
  !
  INTEGER i, j, k, l, N, ii, Nm, jp1
  REAL Ar(Nm,*), Ai(Nm,*), D(*), E(*), E2(*), Tau(2,*)
  REAL f, g, h, fi, gi, hh, si, scale
  REAL PYTHAG
  !
  !***FIRST EXECUTABLE STATEMENT  HTRIDI
  Tau(1,N) = 1.0E0
  Tau(2,N) = 0.0E0
  !
  DO i = 1, N
    D(i) = Ar(i,i)
  ENDDO
  !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO ii = 1, N
    i = N + 1 - ii
    l = i - 1
    h = 0.0E0
    scale = 0.0E0
    IF ( l>=1 ) THEN
      !     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
      DO k = 1, l
        scale = scale + ABS(Ar(i,k)) + ABS(Ai(i,k))
      ENDDO
      !
      IF ( scale/=0.0E0 ) THEN
        !
        DO k = 1, l
          Ar(i,k) = Ar(i,k)/scale
          Ai(i,k) = Ai(i,k)/scale
          h = h + Ar(i,k)*Ar(i,k) + Ai(i,k)*Ai(i,k)
        ENDDO
        !
        E2(i) = scale*scale*h
        g = SQRT(h)
        E(i) = scale*g
        f = PYTHAG(Ar(i,l),Ai(i,l))
        !     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
        IF ( f==0.0E0 ) THEN
          Tau(1,l) = -Tau(1,i)
          si = Tau(2,i)
          Ar(i,l) = g
          GOTO 50
        ELSE
          Tau(1,l) = (Ai(i,l)*Tau(2,i)-Ar(i,l)*Tau(1,i))/f
          si = (Ar(i,l)*Tau(2,i)+Ai(i,l)*Tau(1,i))/f
          h = h + f*g
          g = 1.0E0 + g/f
          Ar(i,l) = g*Ar(i,l)
          Ai(i,l) = g*Ai(i,l)
          IF ( l/=1 ) GOTO 50
          GOTO 100
        ENDIF
      ELSE
        Tau(1,l) = 1.0E0
        Tau(2,l) = 0.0E0
      ENDIF
    ENDIF
    E(i) = 0.0E0
    E2(i) = 0.0E0
    GOTO 150
    50     f = 0.0E0
    !
    DO j = 1, l
      g = 0.0E0
      gi = 0.0E0
      !     .......... FORM ELEMENT OF A*U ..........
      DO k = 1, j
        g = g + Ar(j,k)*Ar(i,k) + Ai(j,k)*Ai(i,k)
        gi = gi - Ar(j,k)*Ai(i,k) + Ai(j,k)*Ar(i,k)
      ENDDO
      !
      jp1 = j + 1
      IF ( l>=jp1 ) THEN
        !
        DO k = jp1, l
          g = g + Ar(k,j)*Ar(i,k) - Ai(k,j)*Ai(i,k)
          gi = gi - Ar(k,j)*Ai(i,k) - Ai(k,j)*Ar(i,k)
        ENDDO
      ENDIF
      !     .......... FORM ELEMENT OF P ..........
      E(j) = g/h
      Tau(2,j) = gi/h
      f = f + E(j)*Ar(i,j) - Tau(2,j)*Ai(i,j)
    ENDDO
    !
    hh = f/(h+h)
    !     .......... FORM REDUCED A ..........
    DO j = 1, l
      f = Ar(i,j)
      g = E(j) - hh*f
      E(j) = g
      fi = -Ai(i,j)
      gi = Tau(2,j) - hh*fi
      Tau(2,j) = -gi
      !
      DO k = 1, j
        Ar(j,k) = Ar(j,k) - f*E(k) - g*Ar(i,k) + fi*Tau(2,k) + gi*Ai(i,k)
        Ai(j,k) = Ai(j,k) - f*Tau(2,k) - g*Ai(i,k) - fi*E(k) - gi*Ar(i,k)
      ENDDO
    ENDDO
    !
    100 CONTINUE
    DO k = 1, l
      Ar(i,k) = scale*Ar(i,k)
      Ai(i,k) = scale*Ai(i,k)
    ENDDO
    !
    Tau(2,l) = -si
    150    hh = D(i)
    D(i) = Ar(i,i)
    Ar(i,i) = hh
    Ai(i,i) = scale*SQRT(h)
  ENDDO
  !
END SUBROUTINE HTRIDI
