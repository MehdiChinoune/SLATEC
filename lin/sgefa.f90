!*==SGEFA.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SGEFA
SUBROUTINE SGEFA(A,Lda,N,Ipvt,Info)
  IMPLICIT NONE
  !*--SGEFA5
  !***BEGIN PROLOGUE  SGEFA
  !***PURPOSE  Factor a matrix using Gaussian elimination.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2A1
  !***TYPE      SINGLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
  !***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
  !             MATRIX FACTORIZATION
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     SGEFA factors a real matrix by Gaussian elimination.
  !
  !     SGEFA is usually called by SGECO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (Time for SGECO) = (1 + 9/N)*(Time for SGEFA) .
  !
  !     On Entry
  !
  !        A       REAL(LDA, N)
  !                the matrix to be factored.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        A       an upper triangular matrix and the multipliers
  !                which were used to obtain it.
  !                The factorization can be written  A = L*U, where
  !                L  is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that SGESL or SGEDI will divide by zero
  !                     if called.  Use  RCOND  in SGECO for a reliable
  !                     indication of singularity.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  ISAMAX, SAXPY, SSCAL
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGEFA
  INTEGER Lda, N, Ipvt(*), Info
  REAL A(Lda,*)
  !
  REAL t
  INTEGER ISAMAX, j, k, kp1, l, nm1
  !
  !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  !
  !***FIRST EXECUTABLE STATEMENT  SGEFA
  Info = 0
  nm1 = N - 1
  IF ( nm1>=1 ) THEN
    DO k = 1, nm1
      kp1 = k + 1
      !
      !        FIND L = PIVOT INDEX
      !
      l = ISAMAX(N-k+1,A(k,k),1) + k - 1
      Ipvt(k) = l
      !
      !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
      !
      IF ( A(l,k)==0.0E0 ) THEN
        Info = k
      ELSE
        !
        !           INTERCHANGE IF NECESSARY
        !
        IF ( l/=k ) THEN
          t = A(l,k)
          A(l,k) = A(k,k)
          A(k,k) = t
        ENDIF
        !
        !           COMPUTE MULTIPLIERS
        !
        t = -1.0E0/A(k,k)
        CALL SSCAL(N-k,t,A(k+1,k),1)
        !
        !           ROW ELIMINATION WITH COLUMN INDEXING
        !
        DO j = kp1, N
          t = A(l,j)
          IF ( l/=k ) THEN
            A(l,j) = A(k,j)
            A(k,j) = t
          ENDIF
          CALL SAXPY(N-k,t,A(k+1,k),1,A(k+1,j),1)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  Ipvt(N) = N
  IF ( A(N,N)==0.0E0 ) Info = N
END SUBROUTINE SGEFA
