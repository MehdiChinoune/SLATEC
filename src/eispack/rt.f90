!** RT
SUBROUTINE RT(Nm,N,A,W,Matz,Z,Fv1,Ierr)
  IMPLICIT NONE
  !>
  !***
  !  Compute the eigenvalues and eigenvectors of a special real
  !            tridiagonal matrix.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A5
  !***
  ! **Type:**      SINGLE PRECISION (RT-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine calls the recommended sequence of subroutines
  !     from the eigensystem subroutine package (EISPACK) to find the
  !     eigenvalues and eigenvectors (if desired) of a special REAL
  !     TRIDIAGONAL matrix.  The property of the matrix required for use
  !     of this subroutine is that the products of pairs of corresponding
  !     off-diagonal elements be all non-negative.  If eigenvectors are
  !     desired, no product can be zero unless both factors are zero.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        A contains the special real tridiagonal matrix in its first
  !          three columns.  The subdiagonal elements are stored in the
  !          last N-1 positions of the first column, the diagonal elements
  !          in the second column, and the superdiagonal elements in the
  !          first N-1 positions of the third column.  Elements A(1,1) and
  !          A(N,3) are arbitrary.  A is a two-dimensional REAL array,
  !          dimensioned A(NM,3).
  !
  !        MATZ is an INTEGER variable set equal to zero if only
  !          eigenvalues are desired.  Otherwise, it is set to any
  !          non-zero integer for both eigenvalues and eigenvectors.
  !
  !     On Output
  !
  !        W contains the eigenvalues in ascending order.  W is a
  !          one-dimensional REAL array, dimensioned W(N).
  !
  !        Z contains the eigenvectors if MATZ is not zero.  The eigen-
  !          vectors are not normalized.  Z is a two-dimensional REAL
  !          array, dimensioned Z(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          10*N       if N is greater than NM,
  !          N+J        if A(J,1)*A(J-1,3) is negative,
  !          2*N+J      if the product is zero with one factor non-zero,
  !                     and MATZ is non-zero;
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !                     The eigenvalues and eigenvectors in the W and Z
  !                     arrays should be correct for indices
  !                     1, 2, ..., IERR-1.
  !
  !        FV1 is a one-dimensional REAL array used for temporary storage,
  !          dimensioned FV1(N).
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
  ! **Routines called:**  FIGI, FIGI2, IMTQL1, IMTQL2

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  !
  INTEGER N, Nm, Ierr, Matz
  REAL A(Nm,3), W(*), Z(Nm,*), Fv1(*)
  !
  !* FIRST EXECUTABLE STATEMENT  RT
  IF ( N>Nm ) THEN
    Ierr = 10*N
    !
  ELSEIF ( Matz/=0 ) THEN
    !     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
    CALL FIGI2(Nm,N,A,W,Fv1,Z,Ierr)
    IF ( Ierr==0 ) CALL IMTQL2(Nm,N,W,Fv1,Z,Ierr)
  ELSE
    !     .......... FIND EIGENVALUES ONLY ..........
    CALL FIGI(Nm,N,A,W,Fv1,Fv1,Ierr)
    IF ( Ierr<=0 ) CALL IMTQL1(N,W,Fv1,Ierr)
  ENDIF
END SUBROUTINE RT
