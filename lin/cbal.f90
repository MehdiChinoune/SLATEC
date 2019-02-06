!*==CBAL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CBAL
      SUBROUTINE CBAL(Nm,N,Ar,Ai,Low,Igh,Scale)
      IMPLICIT NONE
!*--CBAL5
!***BEGIN PROLOGUE  CBAL
!***PURPOSE  Balance a complex general matrix and isolate eigenvalues
!            whenever possible.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1A
!***TYPE      COMPLEX (BALANC-S, CBAL-C)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure
!     CBALANCE, which is a complex version of BALANCE,
!     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
!
!     This subroutine balances a COMPLEX matrix and isolates
!     eigenvalues whenever possible.
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
!        AR and AI contain the real and imaginary parts,
!          respectively, of the complex matrix to be balanced.
!          AR and AI are two-dimensional REAL arrays, dimensioned
!          AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        AR and AI contain the real and imaginary parts,
!          respectively, of the balanced matrix.
!
!        LOW and IGH are two INTEGER variables such that AR(I,J)
!          and AI(I,J) are equal to zero if
!           (1) I is greater than J and
!           (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!        SCALE contains information determining the permutations and
!          scaling factors used.  SCALE is a one-dimensional REAL array,
!          dimensioned SCALE(N).
!
!     Suppose that the principal submatrix in rows LOW through IGH
!     has been balanced, that P(J) denotes the index interchanged
!     with J during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by D(I,J).  Then
!        SCALE(J) = P(J),    for J = 1,...,LOW-1
!                 = D(J,J)       J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     The order in which the interchanges are made is N to IGH+1,
!     then 1 to LOW-1.
!
!     Note that 1 is returned for IGH if IGH is zero formally.
!
!     The ALGOL procedure EXC contained in CBALANCE appears in
!     CBAL  in line.  (Note that the ALGOL roles of identifiers
!     K,L have been reversed.)
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CBAL
!
      INTEGER i , j , k , l , m , N , jj , Nm , Igh , Low , iexc
      REAL Ar(Nm,*) , Ai(Nm,*) , Scale(*)
      REAL c , f , g , r , s , b2 , radix
      LOGICAL noconv
!
!     THE FOLLOWING PORTABLE VALUE OF RADIX WORKS WELL ENOUGH
!     FOR ALL MACHINES WHOSE BASE IS A POWER OF TWO.
!
!***FIRST EXECUTABLE STATEMENT  CBAL
      radix = 16
!
      b2 = radix*radix
      k = 1
      l = N
      GOTO 200
!     .......... IN-LINE PROCEDURE FOR ROW AND
!                COLUMN EXCHANGE ..........
 100  Scale(m) = j
      IF ( j/=m ) THEN
!
        DO i = 1 , l
          f = Ar(i,j)
          Ar(i,j) = Ar(i,m)
          Ar(i,m) = f
          f = Ai(i,j)
          Ai(i,j) = Ai(i,m)
          Ai(i,m) = f
        ENDDO
!
        DO i = k , N
          f = Ar(j,i)
          Ar(j,i) = Ar(m,i)
          Ar(m,i) = f
          f = Ai(j,i)
          Ai(j,i) = Ai(m,i)
          Ai(m,i) = f
        ENDDO
      ENDIF
!
      IF ( iexc==2 ) THEN
!     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT ..........
        k = k + 1
        GOTO 400
      ELSE
!     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                AND PUSH THEM DOWN ..........
        IF ( l==1 ) GOTO 600
        l = l - 1
      ENDIF
!     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
 200  DO jj = 1 , l
        j = l + 1 - jj
!
        DO i = 1 , l
          IF ( i/=j ) THEN
            IF ( Ar(j,i)/=0.0E0.OR.Ai(j,i)/=0.0E0 ) GOTO 300
          ENDIF
        ENDDO
!
        m = l
        iexc = 1
        GOTO 100
!
 300  ENDDO
!
 400  DO j = k , l
!
        DO i = k , l
          IF ( i/=j ) THEN
            IF ( Ar(i,j)/=0.0E0.OR.Ai(i,j)/=0.0E0 ) GOTO 500
          ENDIF
        ENDDO
!
        m = k
        iexc = 2
        GOTO 100
 500  ENDDO
!     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO i = k , l
        Scale(i) = 1.0E0
      ENDDO
      DO
!     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
        noconv = .FALSE.
!
        DO i = k , l
          c = 0.0E0
          r = 0.0E0
!
          DO j = k , l
            IF ( j/=i ) THEN
              c = c + ABS(Ar(j,i)) + ABS(Ai(j,i))
              r = r + ABS(Ar(i,j)) + ABS(Ai(i,j))
            ENDIF
          ENDDO
!     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
          IF ( c/=0.0E0.AND.r/=0.0E0 ) THEN
            g = r/radix
            f = 1.0E0
            s = c + r
            DO WHILE ( c<g )
              f = f*radix
              c = c*b2
            ENDDO
            g = r*radix
            DO WHILE ( c>=g )
              f = f/radix
              c = c/b2
            ENDDO
!     .......... NOW BALANCE ..........
            IF ( (c+r)/f<0.95E0*s ) THEN
              g = 1.0E0/f
              Scale(i) = Scale(i)*f
              noconv = .TRUE.
!
              DO j = k , N
                Ar(i,j) = Ar(i,j)*g
                Ai(i,j) = Ai(i,j)*g
              ENDDO
!
              DO j = 1 , l
                Ar(j,i) = Ar(j,i)*f
                Ai(j,i) = Ai(j,i)*f
              ENDDO
            ENDIF
          ENDIF
!
        ENDDO
!
        IF ( .NOT.(noconv) ) EXIT
      ENDDO
!
 600  Low = k
      Igh = l
      END SUBROUTINE CBAL
