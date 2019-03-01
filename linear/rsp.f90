!DECK RSP
SUBROUTINE RSP(Nm,N,Nv,A,W,Matz,Z,Fv1,Fv2,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  RSP
  !***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a real symmetric matrix packed into a one dimensional
  !            array.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4A1
  !***TYPE      SINGLE PRECISION (RSP-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine calls the recommended sequence of
  !     subroutines from the eigensystem subroutine package (EISPACK)
  !     to find the eigenvalues and eigenvectors (if desired)
  !     of a REAL SYMMETRIC PACKED matrix.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, Z, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        NV is an INTEGER variable set equal to the dimension of the
  !          array A as specified in the calling program.  NV must not
  !          be less than  N*(N+1)/2.
  !
  !        A contains the lower triangle, stored row-wise, of the real
  !          symmetric packed matrix.  A is a one-dimensional REAL
  !          array, dimensioned A(NV).
  !
  !        MATZ is an INTEGER variable set equal to zero if only
  !          eigenvalues are desired.  Otherwise, it is set to any
  !          non-zero integer for both eigenvalues and eigenvectors.
  !
  !     On Output
  !
  !        A has been destroyed.
  !
  !        W contains the eigenvalues in ascending order.  W is a
  !          one-dimensional REAL array, dimensioned W(N).
  !
  !        Z contains the eigenvectors if MATZ is not zero.  The eigen-
  !          vectors are orthonormal.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          10*N       if N is greater than NM,
  !          20*N       if NV is less than N*(N+1)/2,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !                     The eigenvalues and eigenvectors in the W and Z
  !                     arrays should be correct for indices
  !                     1, 2, ..., IERR-1.
  !
  !        FV1 and FV2 are one-dimensional REAL arrays used for temporary
  !          storage, dimensioned FV1(N) and FV2(N).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  TQL2, TQLRAT, TRBAK3, TRED3
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  RSP
  !
  INTEGER i, j, N, Nm, Nv, Ierr, Matz
  REAL A(*), W(*), Z(Nm,*), Fv1(*), Fv2(*)
  !
  !***FIRST EXECUTABLE STATEMENT  RSP
  IF ( N>Nm ) THEN
    Ierr = 10*N
  ELSEIF ( Nv>=(N*(N+1))/2 ) THEN
    !
    CALL TRED3(N,Nv,A,W,Fv1,Fv2)
    IF ( Matz/=0 ) THEN
      !     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      DO i = 1, N
        !
        DO j = 1, N
          Z(j,i) = 0.0E0
        ENDDO
        !
        Z(i,i) = 1.0E0
      ENDDO
      !
      CALL TQL2(Nm,N,W,Fv1,Z,Ierr)
      IF ( Ierr==0 ) CALL TRBAK3(Nm,N,Nv,A,N,Z)
    ELSE
      !     .......... FIND EIGENVALUES ONLY ..........
      CALL TQLRAT(N,W,Fv2,Ierr)
    ENDIF
  ELSE
    Ierr = 20*N
  ENDIF
END SUBROUTINE RSP
