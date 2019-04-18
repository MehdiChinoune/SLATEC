!** SSPEV
SUBROUTINE SSPEV(A,N,E,V,Ldv,Work,Job,Info)
  !>
  !***
  !  Compute the eigenvalues and, optionally, the eigenvectors
  !            of a real symmetric matrix stored in packed form.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A1
  !***
  ! **Type:**      SINGLE PRECISION (SSPEV-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK, PACKED, SYMMETRIC
  !***
  ! **Author:**  Kahaner, D. K., (NBS)
  !           Moler, C. B., (U. of New Mexico)
  !           Stewart, G. W., (U. of Maryland)
  !***
  ! **Description:**
  !
  !     Abstract
  !      SSPEV computes the eigenvalues and, optionally, the eigenvectors
  !      of a real symmetric matrix stored in packed form.
  !
  !     Call Sequence Parameters-
  !       (The values of parameters marked with * (star) will be  changed
  !         by SSPEV.)
  !
  !        A*      REAL(N*(N+1)/2)
  !                real symmetric packed input matrix.  Contains upper
  !                triangle and diagonal of A, by column (elements
  !                11, 12, 22, 13, 23, 33, ...).
  !
  !        N       INTEGER
  !                set by the user to
  !                the order of the matrix A.
  !
  !        E*      REAL(N)
  !                on return from SSPEV, E contains the eigenvalues of A.
  !                See also INFO below.
  !
  !        V*      REAL(LDV,N)
  !                on return from SSPEV, if the user has set JOB
  !                = 0        V is not referenced.
  !                = nonzero  the N eigenvectors of A are stored in the
  !                first N columns of V.  See also INFO below.
  !
  !        LDV     INTEGER
  !                set by the user to
  !                the leading dimension of the array V if JOB is also
  !                set nonzero.  In that case, N must be .LE. LDV.
  !                If JOB is set to zero, LDV is not referenced.
  !
  !        WORK*   REAL(2N)
  !                temporary storage vector.  Contents changed by SSPEV.
  !
  !        JOB     INTEGER
  !                set by the user to
  !                = 0        eigenvalues only to be calculated by SSPEV.
  !                           Neither V nor LDV are referenced.
  !                = nonzero  eigenvalues and vectors to be calculated.
  !                           In this case, A & V must be distinct arrays.
  !                           Also, if LDA .GT. LDV, SSPEV changes all the
  !                           elements of A thru column N.  If LDA < LDV,
  !                           SSPEV changes all the elements of V through
  !                           column N.  If LDA=LDV, only A(I,J) and V(I,
  !                           J) for I,J = 1,...,N are changed by SSPEV.
  !
  !       INFO*   INTEGER
  !               on return from SSPEV, the value of INFO is
  !               = 0 for normal return.
  !               = K if the eigenvalue iteration fails to converge.
  !                   Eigenvalues and vectors 1 through K-1 are correct.
  !
  !
  !     Error Messages-
  !          No. 1   recoverable  N is greater than LDV and JOB is nonzero
  !          No. 2   recoverable  N is less than one
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  IMTQL2, TQLRAT, TRBAK3, TRED3, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800808  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG
  INTEGER Job
  INTEGER i, Info, j, Ldv, m, N
  REAL A(*), E(*), V(Ldv,*), Work(*)
  !* FIRST EXECUTABLE STATEMENT  SSPEV
  IF ( N>Ldv ) CALL XERMSG('SLATEC','SSPEV','N .GT. LDV.',1,1)
  IF ( N>Ldv ) RETURN
  IF ( N<1 ) CALL XERMSG('SLATEC','SSPEV','N .LT. 1',2,1)
  IF ( N<1 ) RETURN
  !
  !       CHECK N=1 CASE
  !
  E(1) = A(1)
  Info = 0
  IF ( N==1 ) RETURN
  !
  IF ( Job/=0 ) THEN
    !
    !     EIGENVALUES AND EIGENVECTORS
    !
    CALL TRED3(N,1,A,E,Work(1),Work(1))
    DO i = 1, N
      DO j = 1, N
        V(i,j) = 0.
      END DO
      V(i,i) = 1.
    END DO
    CALL IMTQL2(Ldv,N,E,Work,V,Info)
    m = N
    IF ( Info/=0 ) m = Info - 1
    CALL TRBAK3(Ldv,N,1,A,m,V)
    RETURN
  END IF
  !
  !     EIGENVALUES ONLY
  !
  CALL TRED3(N,1,A,E,Work(1),Work(N+1))
  CALL TQLRAT(N,E,Work(N+1),Info)
  RETURN
END SUBROUTINE SSPEV
