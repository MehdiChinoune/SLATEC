!DECK CORTH
SUBROUTINE CORTH(Nm,N,Low,Igh,Ar,Ai,Ortr,Orti)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CORTH
  !***PURPOSE  Reduce a complex general matrix to complex upper Hessenberg
  !            form using unitary similarity transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B2
  !***TYPE      COMPLEX (ORTHES-S, CORTH-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of a complex analogue of
  !     the ALGOL procedure ORTHES, NUM. MATH. 12, 349-368(1968)
  !     by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     Given a COMPLEX GENERAL matrix, this subroutine
  !     reduces a submatrix situated in rows and columns
  !     LOW through IGH to upper Hessenberg form by
  !     unitary similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR and AI, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A=(AR,AI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the complex input matrix.  AR and AI are two-dimensional
  !          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
  !
  !     On OUTPUT
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the Hessenberg matrix.  Information about the unitary
  !          transformations used in the reduction is stored in the
  !          remaining triangles under the Hessenberg matrix.
  !
  !        ORTR and ORTI contain further information about the unitary
  !          transformations.  Only elements LOW through IGH are used.
  !          ORTR and ORTI are one-dimensional REAL arrays, dimensioned
  !          ORTR(IGH) and ORTI(IGH).
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
  !***END PROLOGUE  CORTH
  !
  INTEGER i, j, m, N, ii, jj, la, mp, Nm, Igh, kp1, Low
  REAL Ar(Nm,*), Ai(Nm,*), Ortr(*), Orti(*)
  REAL f, g, h, fi, fr, scale
  REAL PYTHAG
  !
  !***FIRST EXECUTABLE STATEMENT  CORTH
  la = Igh - 1
  kp1 = Low + 1
  IF ( la>=kp1 ) THEN
    !
    DO m = kp1, la
      h = 0.0E0
      Ortr(m) = 0.0E0
      Orti(m) = 0.0E0
      scale = 0.0E0
      !     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
      DO i = m, Igh
        scale = scale + ABS(Ar(i,m-1)) + ABS(Ai(i,m-1))
      ENDDO
      !
      IF ( scale/=0.0E0 ) THEN
        mp = m + Igh
        !     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
        DO ii = m, Igh
          i = mp - ii
          Ortr(i) = Ar(i,m-1)/scale
          Orti(i) = Ai(i,m-1)/scale
          h = h + Ortr(i)*Ortr(i) + Orti(i)*Orti(i)
        ENDDO
        !
        g = SQRT(h)
        f = PYTHAG(Ortr(m),Orti(m))
        IF ( f==0.0E0 ) THEN
          !
          Ortr(m) = g
          Ar(m,m-1) = scale
        ELSE
          h = h + f*g
          g = g/f
          Ortr(m) = (1.0E0+g)*Ortr(m)
          Orti(m) = (1.0E0+g)*Orti(m)
        ENDIF
        !     .......... FORM (I-(U*UT)/H) * A ..........
        DO j = m, N
          fr = 0.0E0
          fi = 0.0E0
          !     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
          DO ii = m, Igh
            i = mp - ii
            fr = fr + Ortr(i)*Ar(i,j) + Orti(i)*Ai(i,j)
            fi = fi + Ortr(i)*Ai(i,j) - Orti(i)*Ar(i,j)
          ENDDO
          !
          fr = fr/h
          fi = fi/h
          !
          DO i = m, Igh
            Ar(i,j) = Ar(i,j) - fr*Ortr(i) + fi*Orti(i)
            Ai(i,j) = Ai(i,j) - fr*Orti(i) - fi*Ortr(i)
          ENDDO
          !
        ENDDO
        !     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
        DO i = 1, Igh
          fr = 0.0E0
          fi = 0.0E0
          !     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
          DO jj = m, Igh
            j = mp - jj
            fr = fr + Ortr(j)*Ar(i,j) - Orti(j)*Ai(i,j)
            fi = fi + Ortr(j)*Ai(i,j) + Orti(j)*Ar(i,j)
          ENDDO
          !
          fr = fr/h
          fi = fi/h
          !
          DO j = m, Igh
            Ar(i,j) = Ar(i,j) - fr*Ortr(j) - fi*Orti(j)
            Ai(i,j) = Ai(i,j) + fr*Orti(j) - fi*Ortr(j)
          ENDDO
          !
        ENDDO
        !
        Ortr(m) = scale*Ortr(m)
        Orti(m) = scale*Orti(m)
        Ar(m,m-1) = -g*Ar(m,m-1)
        Ai(m,m-1) = -g*Ai(m,m-1)
      ENDIF
    ENDDO
  ENDIF
  !
END SUBROUTINE CORTH
