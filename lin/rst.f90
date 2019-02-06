!*==RST.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK RST
      SUBROUTINE RST(Nm,N,W,E,Matz,Z,Ierr)
      IMPLICIT NONE
!*--RST5
!***BEGIN PROLOGUE  RST
!***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
!            of a real symmetric tridiagonal matrix.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5
!***TYPE      SINGLE PRECISION (RST-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a REAL SYMMETRIC TRIDIAGONAL matrix.
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
!        W contains the diagonal elements of the real symmetric
!          tridiagonal matrix.  W is a one-dimensional REAL array,
!          dimensioned W(N).
!
!        E contains the subdiagonal elements of the matrix in its last
!          N-1 positions.  E(1) is arbitrary.  E is a one-dimensional
!          REAL array, dimensioned E(N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        W contains the eigenvalues in ascending order.
!
!        Z contains the eigenvectors if MATZ is not zero.  The eigen-
!          vectors are orthonormal.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues and eigenvectors in the W and Z
!                     arrays should be correct for indices
!                     1, 2, ..., IERR-1.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  IMTQL1, IMTQL2
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RST
!
      INTEGER i , j , N , Nm , Ierr , Matz
      REAL W(*) , E(*) , Z(Nm,*)
!
!***FIRST EXECUTABLE STATEMENT  RST
      IF ( N>Nm ) THEN
        Ierr = 10*N
!
      ELSEIF ( Matz/=0 ) THEN
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
        DO i = 1 , N
!
          DO j = 1 , N
            Z(j,i) = 0.0E0
          ENDDO
!
          Z(i,i) = 1.0E0
        ENDDO
!
        CALL IMTQL2(Nm,N,W,E,Z,Ierr)
      ELSE
!     .......... FIND EIGENVALUES ONLY ..........
        CALL IMTQL1(N,W,E,Ierr)
      ENDIF
      END SUBROUTINE RST
