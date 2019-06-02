!** CGEEV
SUBROUTINE CGEEV(A,Lda,N,E,V,Ldv,Work,Job,Info)
  !>
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a complex general matrix.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D4A4
  !***
  ! **Type:**      COMPLEX (SGEEV-S, CGEEV-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
  !***
  ! **Author:**  Kahaner, D. K., (NBS)
  !           Moler, C. B., (U. of New Mexico)
  !           Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     Abstract
  !      CGEEV computes the eigenvalues and, optionally,
  !      the eigenvectors of a general complex matrix.
  !
  !     Call Sequence Parameters-
  !       (The values of parameters marked with * (star) will be changed
  !         by CGEEV.)
  !
  !        A*      COMPLEX(LDA,N)
  !                complex nonsymmetric input matrix.
  !
  !        LDA     INTEGER
  !                set by the user to
  !                the leading dimension of the complex array A.
  !
  !        N       INTEGER
  !                set by the user to
  !                the order of the matrices A and V, and
  !                the number of elements in E.
  !
  !        E*      COMPLEX(N)
  !                on return from CGEEV E contains the eigenvalues of A.
  !                See also INFO below.
  !
  !        V*      COMPLEX(LDV,N)
  !                on return from CGEEV if the user has set JOB
  !                = 0        V is not referenced.
  !                = nonzero  the N eigenvectors of A are stored in the
  !                first N columns of V.  See also INFO below.
  !                (If the input matrix A is nearly degenerate, V
  !                 will be badly conditioned, i.e. have nearly
  !                 dependent columns.)
  !
  !        LDV     INTEGER
  !                set by the user to
  !                the leading dimension of the array V if JOB is also
  !                set nonzero.  In that case N must be .LE. LDV.
  !                If JOB is set to zero LDV is not referenced.
  !
  !        WORK*   REAL(3N)
  !                temporary storage vector.  Contents changed by CGEEV.
  !
  !        JOB     INTEGER
  !                set by the user to
  !                = 0        eigenvalues only to be calculated by CGEEV.
  !                           neither V nor LDV are referenced.
  !                = nonzero  eigenvalues and vectors to be calculated.
  !                           In this case A & V must be distinct arrays.
  !                           Also,  if LDA > LDV,  CGEEV changes all the
  !                           elements of A thru column N.  If LDA < LDV,
  !                           CGEEV changes all the elements of V through
  !                           column N.  If LDA = LDV only A(I,J) and V(I,
  !                           J) for I,J = 1,...,N are changed by CGEEV.
  !
  !        INFO*   INTEGER
  !                on return from CGEEV the value of INFO is
  !                = 0  normal return, calculation successful.
  !                = K  if the eigenvalue iteration fails to converge,
  !                     eigenvalues K+1 through N are correct, but
  !                     no eigenvectors were computed even if they were
  !                     requested (JOB nonzero).
  !
  !      Error Messages
  !           No. 1  recoverable  N is greater than LDA
  !           No. 2  recoverable  N is less than one.
  !           No. 3  recoverable  JOB is nonzero and N is greater than LDV
  !           No. 4  warning      LDA > LDV,  elements of A other than the
  !                               N by N input elements have been changed
  !           No. 5  warning      LDA < LDV,  elements of V other than the
  !                               N by N output elements have been changed
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CBABK2, CBAL, COMQR, COMQR2, CORTH, SCOPY, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG
  INTEGER :: Info, Job, Lda, Ldv, N
  REAL :: A(2*Lda*N), E(2*N), Work(3*N), V(2*Ldv*N)
  INTEGER :: i, ihi, ilo, j, k, l, mdim, m
  !* FIRST EXECUTABLE STATEMENT  CGEEV
  IF ( N>Lda ) CALL XERMSG('CGEEV','N .GT. LDA.',1,1)
  IF ( N>Lda ) RETURN
  IF ( N<1 ) CALL XERMSG('CGEEV','N .LT. 1',2,1)
  IF ( N<1 ) RETURN
  IF ( N/=1.OR.Job/=0 ) THEN
    mdim = 2*Lda
    IF ( Job/=0 ) THEN
      IF ( N>Ldv ) CALL XERMSG('CGEEV','JOB .NE. 0 AND N .GT. LDV.',3,1)
      IF ( N>Ldv ) RETURN
      IF ( N==1 ) GOTO 100
      !
      !       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
      !
      mdim = MIN(mdim,2*Ldv)
      IF ( Lda<Ldv ) CALL XERMSG('CGEEV',&
        'LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.',5,0)
      IF ( Lda>Ldv ) THEN
        CALL XERMSG('CGEEV',&
          'LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ELEMENTS HAVE BEEN CHANGED.',4,0)
        l = N - 1
        DO j = 1, l
          i = 2*N
          m = 1 + j*2*Ldv
          k = 1 + j*2*Lda
          CALL SCOPY(i,A(k),1,A(m),1)
        END DO
      END IF
    END IF
    !
    !     SEPARATE REAL AND IMAGINARY PARTS
    !
    DO j = 1, N
      k = (j-1)*mdim + 1
      l = k + N
      CALL SCOPY(N,A(k+1),2,Work(1),1)
      CALL SCOPY(N,A(k),2,A(k),1)
      CALL SCOPY(N,Work(1),1,A(l),1)
    END DO
    !
    !     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
    !
    CALL CBAL(mdim,N,A(1),A(N+1),ilo,ihi,Work(1))
    CALL CORTH(mdim,N,ilo,ihi,A(1),A(N+1),Work(N+1),Work(2*N+1))
    IF ( Job/=0 ) THEN
      !
      !     EIGENVALUES AND EIGENVECTORS.
      !
      CALL COMQR2(mdim,N,ilo,ihi,Work(N+1),Work(2*N+1),A(1),A(N+1),E(1),&
        E(N+1),V(1),V(N+1),Info)
      IF ( Info==0 ) THEN
        CALL CBABK2(mdim,N,ilo,ihi,Work(1),N,V(1),V(N+1))
        !
        !     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
        !
        DO j = 1, N
          k = (j-1)*mdim + 1
          i = (j-1)*2*Ldv + 1
          l = k + N
          CALL SCOPY(N,V(k),1,Work(1),1)
          CALL SCOPY(N,V(l),1,V(i+1),2)
          CALL SCOPY(N,Work(1),1,V(i),2)
        END DO
      END IF
    ELSE
      !
      !     EIGENVALUES ONLY
      !
      CALL COMQR(mdim,N,ilo,ihi,A(1),A(N+1),E(1),E(N+1),Info)
    END IF
    !
    !     CONVERT EIGENVALUES TO COMPLEX STORAGE.
    !
    CALL SCOPY(N,E(1),1,Work(1),1)
    CALL SCOPY(N,E(N+1),1,E(2),2)
    CALL SCOPY(N,Work(1),1,E(1),2)
    RETURN
  END IF
  !
  !     TAKE CARE OF N=1 CASE
  !
  100  E(1) = A(1)
  E(2) = A(2)
  Info = 0
  IF ( Job==0 ) RETURN
  V(1) = A(1)
  V(2) = A(2)
END SUBROUTINE CGEEV
