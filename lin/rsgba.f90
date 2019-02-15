!DECK RSGBA
SUBROUTINE RSGBA(Nm,N,A,B,W,Matz,Z,Fv1,Fv2,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  RSGBA
  !***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a symmetric generalized eigenproblem.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4B1
  !***TYPE      SINGLE PRECISION (RSGBA-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine calls the recommended sequence of
  !     subroutines from the eigensystem subroutine package (EISPACK)
  !     to find the eigenvalues and eigenvectors (if desired)
  !     for the REAL SYMMETRIC generalized eigenproblem  BAx = (LAMBDA)x.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A, B, and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrices A and B.  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        A contains a real symmetric matrix.  A is a two-dimensional
  !          REAL array, dimensioned A(NM,N).
  !
  !        B contains a positive definite real symmetric matrix.  B is a
  !          two-dimensional REAL array, dimensioned B(NM,N).
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
  !        Z contains the eigenvectors if MATZ is not zero.  Z is a
  !          two-dimensional REAL array, dimensioned Z(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          10*N       if N is greater than NM,
  !          7*N+1      if B is not positive definite,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !                     The eigenvalues should be correct for indices
  !                     1, 2, ..., IERR-1, but no eigenvectors are
  !                     computed.
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
  !***ROUTINES CALLED  REBAKB, REDUC2, TQL2, TQLRAT, TRED1, TRED2
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  RSGBA
  !
  INTEGER N, Nm, Ierr, Matz
  REAL A(Nm,*), B(Nm,*), W(*), Z(Nm,*), Fv1(*), Fv2(*)
  !
  !***FIRST EXECUTABLE STATEMENT  RSGBA
  IF ( N<=Nm ) THEN
    !
    CALL REDUC2(Nm,N,A,B,Fv2,Ierr)
    IF ( Ierr==0 ) THEN
      IF ( Matz/=0 ) THEN
        !     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
        CALL TRED2(Nm,N,A,W,Fv1,Z)
        CALL TQL2(Nm,N,W,Fv1,Z,Ierr)
        IF ( Ierr==0 ) CALL REBAKB(Nm,N,B,Fv2,N,Z)
      ELSE
        !     .......... FIND EIGENVALUES ONLY ..........
        CALL TRED1(Nm,N,A,W,Fv1,Fv2)
        CALL TQLRAT(N,W,Fv2,Ierr)
      ENDIF
    ENDIF
  ELSE
    Ierr = 10*N
  ENDIF
END SUBROUTINE RSGBA
