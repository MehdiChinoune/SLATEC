!*==SGEEV.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SGEEV
SUBROUTINE SGEEV(A,Lda,N,E,V,Ldv,Work,Job,Info)
  IMPLICIT NONE
  !*--SGEEV5
  !*** Start of declarations inserted by SPAG
  INTEGER m
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  SGEEV
  !***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a real general matrix.
  !***LIBRARY   SLATEC
  !***CATEGORY  D4A2
  !***TYPE      SINGLE PRECISION (SGEEV-S, CGEEV-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX
  !***AUTHOR  Kahaner, D. K., (NBS)
  !           Moler, C. B., (U. of New Mexico)
  !           Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     Abstract
  !      SGEEV computes the eigenvalues and, optionally,
  !      the eigenvectors of a general real matrix.
  !
  !     Call Sequence Parameters-
  !       (The values of parameters marked with * (star) will be changed
  !         by SGEEV.)
  !
  !        A*      REAL(LDA,N)
  !                real nonsymmetric input matrix.
  !
  !        LDA     INTEGER
  !                set by the user to
  !                the leading dimension of the real array A.
  !
  !        N       INTEGER
  !                set by the user to
  !                the order of the matrices A and V, and
  !                the number of elements in E.
  !
  !        E*      COMPLEX(N)
  !                on return from SGEEV, E contains the eigenvalues of A.
  !                See also INFO below.
  !
  !        V*      COMPLEX(LDV,N)
  !                on return from SGEEV, if the user has set JOB
  !                = 0        V is not referenced.
  !                = nonzero  the N eigenvectors of A are stored in the
  !                first N columns of V.  See also INFO below.
  !                (Note that if the input matrix A is nearly degenerate,
  !                 V may be badly conditioned, i.e., may have nearly
  !                 dependent columns.)
  !
  !        LDV     INTEGER
  !                set by the user to
  !                the leading dimension of the array V if JOB is also
  !                set nonzero.  In that case, N must be .LE. LDV.
  !                If JOB is set to zero, LDV is not referenced.
  !
  !        WORK*   REAL(2N)
  !                temporary storage vector.  Contents changed by SGEEV.
  !
  !        JOB     INTEGER
  !                set by the user to
  !                = 0        eigenvalues only to be calculated by SGEEV.
  !                           Neither V nor LDV is referenced.
  !                = nonzero  eigenvalues and vectors to be calculated.
  !                           In this case, A & V must be distinct arrays.
  !                           Also, if LDA .GT. LDV, SGEEV changes all the
  !                           elements of A thru column N.  If LDA < LDV,
  !                           SGEEV changes all the elements of V through
  !                           column N. If LDA = LDV, only A(I,J) and V(I,
  !                           J) for I,J = 1,...,N are changed by SGEEV.
  !
  !        INFO*   INTEGER
  !                on return from SGEEV the value of INFO is
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
  !           No. 4  warning      LDA > LDV, elements of A other than the
  !                               N by N input elements have been changed.
  !           No. 5  warning      LDA < LDV, elements of V other than the
  !                               N x N output elements have been changed.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  BALANC, BALBAK, HQR, HQR2, ORTHES, ORTRAN, SCOPY,
  !                    SCOPYM, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  SGEEV
  INTEGER i, ihi, ilo, Info, j, jb, Job, k, km, kp, l, Lda, &
    Ldv, mdim, N
  REAL A(*), E(*), Work(*), V(*)
  !***FIRST EXECUTABLE STATEMENT  SGEEV
  IF ( N>Lda ) CALL XERMSG('SLATEC','SGEEV','N .GT. LDA.',1,1)
  IF ( N>Lda ) RETURN
  IF ( N<1 ) CALL XERMSG('SLATEC','SGEEV','N .LT. 1',2,1)
  IF ( N<1 ) RETURN
  IF ( N/=1.OR.Job/=0 ) THEN
    mdim = Lda
    IF ( Job/=0 ) THEN
      IF ( N>Ldv ) CALL XERMSG('SLATEC','SGEEV','JOB .NE. 0 AND N .GT. LDV.'&
        ,3,1)
      IF ( N>Ldv ) RETURN
      IF ( N==1 ) GOTO 100
      !
      !       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
      !
      mdim = MIN(Lda,Ldv)
      IF ( Lda<Ldv ) CALL XERMSG('SLATEC','SGEEV',&
        'LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT '&
        //'ELEMENTS HAVE BEEN CHANGED.',5,0)
      IF ( Lda>Ldv ) THEN
        CALL XERMSG('SLATEC','SGEEV',&
          'LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT '&
          //'ELEMENTS HAVE BEEN CHANGED.',4,0)
        l = N - 1
        DO j = 1, l
          m = 1 + j*Ldv
          k = 1 + j*Lda
          CALL SCOPY(N,A(k),1,A(m),1)
        ENDDO
      ENDIF
    ENDIF
    !
    !     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
    !
    CALL BALANC(mdim,N,A,ilo,ihi,Work(1))
    CALL ORTHES(mdim,N,ilo,ihi,A,Work(N+1))
    IF ( Job/=0 ) THEN
      !
      !     EIGENVALUES AND EIGENVECTORS.
      !
      CALL ORTRAN(mdim,N,ilo,ihi,A,Work(N+1),V)
      CALL HQR2(mdim,N,ilo,ihi,A,E(1),E(N+1),V,Info)
      IF ( Info==0 ) THEN
        CALL BALBAK(mdim,N,ilo,ihi,Work(1),N,V)
        !
        !     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
        !
        DO jb = 1, N
          j = N + 1 - jb
          i = N + j
          k = (j-1)*mdim + 1
          kp = k + mdim
          km = k - mdim
          IF ( E(i)>=0.0E0 ) CALL SCOPY(N,V(k),1,Work(1),2)
          IF ( E(i)<0.0E0 ) CALL SCOPY(N,V(km),1,Work(1),2)
          IF ( E(i)==0.0E0 ) CALL SCOPY(N,0.0E0,0,Work(2),2)
          IF ( E(i)>0.0E0 ) CALL SCOPY(N,V(kp),1,Work(2),2)
          IF ( E(i)<0.0E0 ) CALL SCOPYM(N,V(k),1,Work(2),2)
          l = 2*(j-1)*Ldv + 1
          CALL SCOPY(2*N,Work(1),1,V(l),1)
        ENDDO
      ENDIF
    ELSE
      !
      !     EIGENVALUES ONLY
      !
      CALL HQR(Lda,N,ilo,ihi,A,E(1),E(N+1),Info)
    ENDIF
    !
    !     CONVERT EIGENVALUES TO COMPLEX STORAGE.
    !
    CALL SCOPY(N,E(1),1,Work(1),1)
    CALL SCOPY(N,E(N+1),1,E(2),2)
    CALL SCOPY(N,Work(1),1,E(1),2)
    RETURN
  ENDIF
  !
  !     TAKE CARE OF N=1 CASE
  !
  100  E(1) = A(1)
  E(2) = 0.E0
  Info = 0
  IF ( Job==0 ) RETURN
  V(1) = A(1)
  V(2) = 0.E0
END SUBROUTINE SGEEV
