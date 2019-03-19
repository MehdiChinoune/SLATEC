!** SSIEV
SUBROUTINE SSIEV(A,Lda,N,E,Work,Job,Info)
  IMPLICIT NONE
  !>
  !***
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a real symmetric matrix.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D4A1
  !***
  ! **Type:**      SINGLE PRECISION (SSIEV-S, CHIEV-C)
  !***
  ! **Keywords:**  COMPLEX HERMITIAN, EIGENVALUES, EIGENVECTORS, MATRIX,
  !             SYMMETRIC
  !***
  ! **Author:**  Kahaner, D. K., (NBS)
  !           Moler, C. B., (U. of New Mexico)
  !           Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     Abstract
  !      SSIEV computes the eigenvalues and, optionally, the eigenvectors
  !      of a real symmetric matrix.
  !
  !     Call Sequence Parameters-
  !       (The values of parameters marked with * (star) will be  changed
  !         by SSIEV.)
  !
  !       A*      REAL (LDA,N)
  !               real symmetric input matrix.
  !               Only the diagonal and upper triangle of A must be input,
  !               as SSIEV copies the upper triangle to the lower.
  !               That is, the user must define A(I,J), I=1,..N, and J=I,.
  !               ..,N.
  !               On return from SSIEV, if the user has set JOB
  !               = 0        the lower triangle of A has been altered.
  !               = nonzero  the N eigenvectors of A are stored in its
  !               first N columns.  See also INFO below.
  !
  !       LDA     INTEGER
  !               set by the user to
  !               the leading dimension of the array A.
  !
  !       N       INTEGER
  !               set by the user to
  !               the order of the matrix A and
  !               the number of elements in E.
  !
  !       E*      REAL (N)
  !               on return from SSIEV, E contains the N
  !               eigenvalues of A.  See also INFO below.
  !
  !       WORK*   REAL (2*N)
  !               temporary storage vector.  Contents changed by SSIEV.
  !
  !       JOB     INTEGER
  !               set by user on input
  !               = 0         only calculate eigenvalues of A.
  !               = nonzero   calculate eigenvalues and eigenvectors of A.
  !
  !       INFO*   INTEGER
  !               on return from SSIEV, the value of INFO is
  !               = 0 for normal return.
  !               = K if the eigenvalue iteration fails to converge.
  !                   eigenvalues and vectors 1 through K-1 are correct.
  !
  !
  !     Error Messages-
  !          No. 1   recoverable  N is greater than LDA
  !          No. 2   recoverable  N is less than one
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  IMTQL2, TQLRAT, TRED1, TRED2, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  
  INTEGER i, j
  INTEGER Info, Job, Lda, N
  REAL A(Lda,*), E(*), Work(*)
  !* FIRST EXECUTABLE STATEMENT  SSIEV
  IF ( N>Lda ) CALL XERMSG('SLATEC','SSIEV','N .GT. LDA.',1,1)
  IF ( N>Lda ) RETURN
  IF ( N<1 ) CALL XERMSG('SLATEC','SSIEV','N .LT. 1',2,1)
  IF ( N<1 ) RETURN
  !
  !       CHECK N=1 CASE
  !
  E(1) = A(1,1)
  Info = 0
  IF ( N==1 ) RETURN
  !
  !     COPY UPPER TRIANGLE TO LOWER
  !
  DO j = 1, N
    DO i = 1, j
      A(j,i) = A(i,j)
    ENDDO
  ENDDO
  !
  IF ( Job/=0 ) THEN
    !
    !     EIGENVALUES AND EIGENVECTORS
    !
    CALL TRED2(Lda,N,A,E,Work,A)
    CALL IMTQL2(Lda,N,E,Work,A,Info)
    RETURN
  ENDIF
  !
  !     EIGENVALUES ONLY
  !
  CALL TRED1(Lda,N,A,E,Work(1),Work(N+1))
  CALL TQLRAT(N,E,Work(N+1),Info)
  RETURN
END SUBROUTINE SSIEV
