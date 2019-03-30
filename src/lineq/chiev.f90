!** CHIEV
SUBROUTINE CHIEV(A,Lda,N,E,V,Ldv,Work,Job,Info)
  IMPLICIT NONE
  !>
  !***
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a complex Hermitian matrix.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D4A3
  !***
  ! **Type:**      COMPLEX (SSIEV-S, CHIEV-C)
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
  !     David Kahaner, Cleve Moler, G. W. Stewart,
  !       N.B.S.         U.N.M.      N.B.S./U.MD.
  !
  !     Abstract
  !      CHIEV computes the eigenvalues and, optionally,
  !      the eigenvectors of a complex Hermitian matrix.
  !
  !     Call Sequence Parameters-
  !       (the values of parameters marked with * (star) will be changed
  !         by CHIEV.)
  !
  !        A*      COMPLEX(LDA,N)
  !                complex Hermitian input matrix.
  !                Only the upper triangle of A need be
  !                filled in.  Elements on diagonal must be real.
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
  !        E*      REAL(N)
  !                on return from CHIEV E contains the eigenvalues of A.
  !                See also INFO below.
  !
  !        V*      COMPLEX(LDV,N)
  !                on return from CHIEV if the user has set JOB
  !                = 0        V is not referenced.
  !                = nonzero  the N eigenvectors of A are stored in the
  !                first N columns of V.  See also INFO below.
  !
  !        LDV     INTEGER
  !                set by the user to
  !                the leading dimension of the array V if JOB is also
  !                set nonzero.  In that case N must be .LE. LDV.
  !                If JOB is set to zero LDV is not referenced.
  !
  !        WORK*   REAL(4N)
  !                temporary storage vector.  Contents changed by CHIEV.
  !
  !        JOB     INTEGER
  !                set by the user to
  !                = 0        eigenvalues only to be calculated by CHIEV.
  !                           Neither V nor LDV are referenced.
  !                = nonzero  eigenvalues and vectors to be calculated.
  !                           In this case A and V must be distinct arrays
  !                           also if LDA .GT. LDV CHIEV changes all the
  !                           elements of A thru column N.  If LDA < LDV
  !                           CHIEV changes all the elements of V through
  !                           column N.  If LDA = LDV only A(I,J) and V(I,
  !                           J) for I,J = 1,...,N are changed by CHIEV.
  !
  !        INFO*   INTEGER
  !                on return from CHIEV the value of INFO is
  !                = 0  normal return, calculation successful.
  !                = K  if the eigenvalue iteration fails to converge,
  !                     eigenvalues (and eigenvectors if requested)
  !                     1 through K-1 are correct.
  !
  !      Error Messages
  !           No. 1  recoverable  N is greater than LDA
  !           No. 2  recoverable  N is less than one.
  !           No. 3  recoverable  JOB is nonzero and N is greater than LDV
  !           No. 4  warning      LDA > LDV,  elements of A other than the
  !                               N by N input elements have been changed
  !           No. 5  warning      LDA < LDV,  elements of V other than the
  !                               N by N output elements have been changed
  !           No. 6  recoverable  nonreal element on diagonal of A.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  HTRIBK, HTRIDI, IMTQL2, SCOPY, SCOPYM, TQLRAT,
  !                    XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  INTEGER i, Info, j, Job, k, l, Lda, Ldv, m, mdim, N
  REAL A(*), E(*), Work(*), V(*)
  !* FIRST EXECUTABLE STATEMENT  CHIEV
  IF ( N>Lda ) CALL XERMSG('SLATEC','CHIEV','N .GT. LDA.',1,1)
  IF ( N>Lda ) RETURN
  IF ( N<1 ) CALL XERMSG('SLATEC','CHIEV','N .LT. 1',2,1)
  IF ( N<1 ) RETURN
  IF ( N/=1.OR.Job/=0 ) THEN
    mdim = 2*Lda
    IF ( Job/=0 ) THEN
      IF ( N>Ldv ) CALL XERMSG('SLATEC','CHIEV','JOB .NE. 0 AND N .GT. LDV.',3,1)
      IF ( N>Ldv ) RETURN
      IF ( N==1 ) GOTO 100
      !
      !       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
      !
      mdim = MIN(mdim,2*Ldv)
      IF ( Lda<Ldv ) CALL XERMSG('SLATEC','CHIEV',&
        'LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.',5,0)
      IF ( Lda>Ldv ) THEN
        CALL XERMSG('SLATEC','CHIEV',&
          'LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ELEMENTS HAVE BEEN CHANGED.',4,0)
        l = N - 1
        DO j = 1, l
          m = 1 + j*2*Ldv
          k = 1 + j*2*Lda
          CALL SCOPY(2*N,A(k),1,A(m),1)
        ENDDO
      ENDIF
    ENDIF
    !
    !     FILL IN LOWER TRIANGLE OF A, COLUMN BY COLUMN.
    !
    DO j = 1, N
      k = (j-1)*(mdim+2) + 1
      IF ( A(k+1)/=0.0 )&
        CALL XERMSG('SLATEC','CHIEV','NONREAL ELEMENT ON DIAGONAL OF A',6,1)
      IF ( A(k+1)/=0.0 ) RETURN
      CALL SCOPY(N-j+1,A(k),mdim,A(k),2)
      CALL SCOPYM(N-j+1,A(k+1),mdim,A(k+1),2)
    ENDDO
    !
    !     SEPARATE REAL AND IMAGINARY PARTS
    !
    DO j = 1, N
      k = (j-1)*mdim + 1
      l = k + N
      CALL SCOPY(N,A(k+1),2,Work(1),1)
      CALL SCOPY(N,A(k),2,A(k),1)
      CALL SCOPY(N,Work(1),1,A(l),1)
    ENDDO
    !
    !    REDUCE A TO TRIDIAGONAL MATRIX.
    !
    CALL HTRIDI(mdim,N,A(1),A(N+1),E,Work(1),Work(N+1),Work(2*N+1))
    IF ( Job/=0 ) THEN
      !
      !     EIGENVALUES AND EIGENVECTORS.
      !
      DO j = 1, N
        k = (j-1)*mdim + 1
        m = k + N - 1
        DO i = k, m
          V(i) = 0.
        ENDDO
        i = k + j - 1
        V(i) = 1.
      ENDDO
      CALL IMTQL2(mdim,N,E,Work(1),V,Info)
      IF ( Info/=0 ) RETURN
      CALL HTRIBK(mdim,N,A(1),A(N+1),Work(2*N+1),N,V(1),V(N+1))
      !
      !    CONVERT EIGENVECTORS TO COMPLEX STORAGE.
      !
      DO j = 1, N
        k = (j-1)*mdim + 1
        i = (j-1)*2*Ldv + 1
        l = k + N
        CALL SCOPY(N,V(k),1,Work(1),1)
        CALL SCOPY(N,V(l),1,V(i+1),2)
        CALL SCOPY(N,Work(1),1,V(i),2)
      ENDDO
      RETURN
    ELSE
      !
      !     EIGENVALUES ONLY.
      !
      CALL TQLRAT(N,E,Work(N+1),Info)
      RETURN
    ENDIF
  ENDIF
  !
  !     TAKE CARE OF N=1 CASE.
  !
  100 CONTINUE
  IF ( A(2)/=0. ) CALL XERMSG('SLATEC','CHIEV',&
    'NONREAL ELEMENT ON DIAGONAL OF A',6,1)
  IF ( A(2)/=0. ) RETURN
  E(1) = A(1)
  Info = 0
  IF ( Job==0 ) RETURN
  V(1) = A(1)
  V(2) = 0.
END SUBROUTINE CHIEV
