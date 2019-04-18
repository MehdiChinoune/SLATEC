!** CH
SUBROUTINE CH(Nm,N,Ar,Ai,W,Matz,Zr,Zi,Fv1,Fv2,Fm1,Ierr)
  !>
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a complex Hermitian matrix.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A3
  !***
  ! **Type:**      COMPLEX (RS-S, CH-C)
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
  !     of a COMPLEX HERMITIAN matrix.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR, AI, ZR and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        N is the order of the matrix A=(AR,AI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the complex Hermitian matrix.  AR and AI are
  !          two-dimensional REAL arrays, dimensioned AR(NM,N)
  !          and AI(NM,N).
  !
  !        MATZ is an INTEGER variable set equal to zero if only
  !          eigenvalues are desired.  Otherwise, it is set to any
  !          non-zero integer for both eigenvalues and eigenvectors.
  !
  !     On OUTPUT
  !
  !        W contains the eigenvalues in ascending order.
  !          W is a one-dimensional REAL array, dimensioned W(N).
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors if MATZ is not zero.  ZR and ZI are
  !          two-dimensional REAL arrays, dimensioned ZR(NM,N) and
  !          ZI(NM,N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          10*N       if N is greater than NM,
  !          J          if the J-th eigenvalue has not been
  !                     determined after a total of 30 iterations.
  !                     The eigenvalues should be correct for indices
  !                     1, 2, ..., IERR-1, but no eigenvectors are
  !                     computed.
  !
  !        FV1 and FV2 are one-dimensional REAL arrays used for
  !          temporary storage, dimensioned FV1(N) and FV2(N).
  !
  !        FM1 is a two-dimensional REAL array used for temporary
  !          storage, dimensioned FM1(2,N).
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
  ! **Routines called:**  HTRIBK, HTRIDI, TQL2, TQLRAT

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, N, Nm, Ierr, Matz
  REAL Ar(Nm,*), Ai(Nm,*), W(*), Zr(Nm,*), Zi(Nm,*)
  REAL Fv1(*), Fv2(*), Fm1(2,*)
  !
  !* FIRST EXECUTABLE STATEMENT  CH
  IF ( N<=Nm ) THEN
    !
    CALL HTRIDI(Nm,N,Ar,Ai,W,Fv1,Fv2,Fm1)
    IF ( Matz/=0 ) THEN
      !     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      DO i = 1, N
        !
        DO j = 1, N
          Zr(j,i) = 0.0E0
        END DO
        !
        Zr(i,i) = 1.0E0
      END DO
      !
      CALL TQL2(Nm,N,W,Fv1,Zr,Ierr)
      IF ( Ierr==0 ) CALL HTRIBK(Nm,N,Ar,Ai,Fm1,N,Zr,Zi)
    ELSE
      !     .......... FIND EIGENVALUES ONLY ..........
      CALL TQLRAT(N,W,Fv2,Ierr)
    END IF
  ELSE
    Ierr = 10*N
  END IF
END SUBROUTINE CH
