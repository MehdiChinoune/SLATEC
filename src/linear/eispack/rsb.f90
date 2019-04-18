!** RSB
SUBROUTINE RSB(Nm,N,Mb,A,W,Matz,Z,Fv1,Fv2,Ierr)
  !>
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a symmetric band matrix.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A6
  !***
  ! **Type:**      SINGLE PRECISION (RSB-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine calls the recommended sequence of
  !     subroutines from the eigensystem subroutine package (EISPACK)
  !     to find the eigenvalues and eigenvectors (if desired)
  !     of a REAL SYMMETRIC BAND matrix.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        MB is the half band width of the matrix, defined as the
  !          number of adjacent diagonals, including the principal
  !          diagonal, required to specify the non-zero portion of the
  !          lower triangle of the matrix.  MB must be less than or
  !          equal to N.  MB is an INTEGER variable.
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
  !        MATZ is an INTEGER variable set equal to zero if only
  !          eigenvalues are desired.  Otherwise, it is set to any
  !          non-zero integer for both eigenvalues and eigenvectors.
  !
  !     On Output
  !
  !        A has been destroyed.
  !
  !        W contains the eigenvalues in ascending order.  W is a one-
  !          dimensional REAL array, dimensioned W(N).
  !
  !        Z contains the eigenvectors if MATZ is not zero.  The
  !          eigenvectors are orthonormal.  Z is a two-dimensional
  !          REAL array, dimensioned Z(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          10*N       if N is greater than NM,
  !          12*N       if MB is either non-positive or greater than N,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !                     The eigenvalues and eigenvectors, if requested,
  !                     should be correct for indices 1, 2, ..., IERR-1.
  !
  !        FV1 and FV2 are one-dimensional REAL arrays used for temporary
  !          storage, dimensioned FV1(N) and FV2(N).
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
  ! **Routines called:**  BANDR, TQL2, TQLRAT

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER N, Mb, Nm, Ierr, Matz
  REAL A(Nm,*), W(*), Z(Nm,*), Fv1(*), Fv2(*)
  LOGICAL tf
  !
  !* FIRST EXECUTABLE STATEMENT  RSB
  IF ( N>Nm ) THEN
    Ierr = 10*N
  ELSEIF ( Mb<=0 ) THEN
    Ierr = 12*N
  ELSEIF ( Mb>N ) THEN
    Ierr = 12*N
    !
  ELSEIF ( Matz/=0 ) THEN
    !     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
    tf = .TRUE.
    CALL BANDR(Nm,N,Mb,A,W,Fv1,Fv1,tf,Z)
    CALL TQL2(Nm,N,W,Fv1,Z,Ierr)
  ELSE
    !     .......... FIND EIGENVALUES ONLY ..........
    tf = .FALSE.
    CALL BANDR(Nm,N,Mb,A,W,Fv1,Fv2,tf,Z)
    CALL TQLRAT(N,W,Fv2,Ierr)
  END IF
END SUBROUTINE RSB
