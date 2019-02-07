!*==ORTHOR.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK ORTHOR
SUBROUTINE ORTHOR(A,N,M,Nrda,Iflag,Irank,Iscale,Diag,Kpivot,Scales,Rows,&
    Rs)
  IMPLICIT NONE
  !*--ORTHOR6
  !*** Start of declarations inserted by SPAG
  REAL A, acc, akk, anorm, as, asave, Diag, diagk, dum, R1MACH, &
    Rows, Rs, rss, sad, Scales, SDOT, sig, sigma, sruro, uro
  INTEGER Iflag, Irank, Iscale, j, jrow, k, kp, Kpivot, l, M, mk, &
    N, Nrda
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  ORTHOR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (ORTHOR-S, DORTHR-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   Reduction of the matrix A to lower triangular form by a sequence of
  !   orthogonal HOUSEHOLDER transformations post-multiplying A
  !
  !   Modeled after the ALGOL codes in the articles in the REFERENCES
  !   section.
  !
  ! **********************************************************************
  !   INPUT
  ! **********************************************************************
  !
  !     A -- Contains the matrix to be decomposed, must be dimensioned
  !           NRDA by N
  !     N -- Number of rows in the matrix, N greater or equal to 1
  !     M -- Number of columns in the matrix, M greater or equal to N
  !     IFLAG -- Indicates the uncertainty in the matrix data
  !             = 0 when the data is to be treated as exact
  !             =-K when the data is assumed to be accurate to about
  !                 K digits
  !     ISCALE -- Scaling indicator
  !               =-1 if the matrix is to be pre-scaled by
  !               columns when appropriate.
  !               Otherwise no scaling will be attempted
  !     NRDA -- Row dimension of A, NRDA greater or equal to N
  !     DIAG,KPIVOT,ROWS -- Arrays of length at least N used internally
  !         ,RS,SCALES         (except for SCALES which is M)
  !
  ! **********************************************************************
  !   OUTPUT
  ! **********************************************************************
  !
  !     IFLAG - status indicator
  !            =1 for successful decomposition
  !            =2 if improper input is detected
  !            =3 if rank of the matrix is less than N
  !     A -- contains the reduced matrix in the strictly lower triangular
  !          part and transformation information
  !     IRANK -- contains the numerically determined matrix rank
  !     DIAG -- contains the diagonal elements of the reduced
  !             triangular matrix
  !     KPIVOT -- Contains the pivotal information, the column
  !               interchanges performed on the original matrix are
  !               recorded here.
  !     SCALES -- contains the column scaling parameters
  !
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***REFERENCES  G. Golub, Numerical methods for solving linear least
  !                 squares problems, Numerische Mathematik 7, (1965),
  !                 pp. 206-216.
  !               P. Businger and G. Golub, Linear least squares
  !                 solutions by Householder transformations, Numerische
  !                 Mathematik  7, (1965), pp. 269-276.
  !***ROUTINES CALLED  CSCALE, R1MACH, SDOT, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  ORTHOR
  DIMENSION A(Nrda,*), Diag(*), Kpivot(*), Rows(*), Rs(*), Scales(*)
  !
  ! END OF ABSTRACT
  !
  ! **********************************************************************
  !
  !     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
  !     BY THE FUNCTION R1MACH.
  !
  ! **********************************************************************
  !
  !***FIRST EXECUTABLE STATEMENT  ORTHOR
  uro = R1MACH(4)
  IF ( M>=N.AND.N>=1.AND.Nrda>=N ) THEN
    !
    acc = 10.*uro
    IF ( Iflag<0 ) acc = MAX(acc,10.**Iflag)
    sruro = SQRT(uro)
    Iflag = 1
    Irank = N
    !
    !     COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM
    !
    anorm = 0.
    DO j = 1, N
      Kpivot(j) = j
      Rows(j) = SDOT(M,A(j,1),Nrda,A(j,1),Nrda)
      Rs(j) = Rows(j)
      anorm = anorm + Rows(j)
    ENDDO
    !
    !     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
    !
    CALL CSCALE(A,Nrda,N,M,Scales,dum,Rows,Rs,anorm,Scales,Iscale,1)
    !
    anorm = SQRT(anorm)
    !
    !
    !     CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF
    !     ORTHOGONAL TRANSFORMATIONS
    !
    !
    DO k = 1, N
      mk = M - k + 1
      IF ( k/=N ) THEN
        kp = k + 1
        !
        !        SEARCHING FOR PIVOTAL ROW
        !
        DO j = k, N
          IF ( Rows(j)<sruro*Rs(j) ) THEN
            Rows(j) = SDOT(mk,A(j,k),Nrda,A(j,k),Nrda)
            Rs(j) = Rows(j)
          ENDIF
          IF ( j/=k ) THEN
            IF ( sigma>=0.99*Rows(j) ) CYCLE
          ENDIF
          sigma = Rows(j)
          jrow = j
        ENDDO
        IF ( jrow/=k ) THEN
          !
          !        PERFORM ROW INTERCHANGE
          !
          l = Kpivot(k)
          Kpivot(k) = Kpivot(jrow)
          Kpivot(jrow) = l
          Rows(jrow) = Rows(k)
          Rows(k) = sigma
          rss = Rs(k)
          Rs(k) = Rs(jrow)
          Rs(jrow) = rss
          DO l = 1, M
            asave = A(k,l)
            A(k,l) = A(jrow,l)
            A(jrow,l) = asave
          ENDDO
        ENDIF
      ENDIF
      !
      !        CHECK RANK OF THE MATRIX
      !
      sig = SDOT(mk,A(k,k),Nrda,A(k,k),Nrda)
      diagk = SQRT(sig)
      IF ( diagk>acc*anorm ) THEN
        !
        !        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
        !
        akk = A(k,k)
        IF ( akk>0. ) diagk = -diagk
        Diag(k) = diagk
        A(k,k) = akk - diagk
        IF ( k/=N ) THEN
          sad = diagk*akk - sig
          DO j = kp, N
            as = SDOT(mk,A(k,k),Nrda,A(j,k),Nrda)/sad
            DO l = k, M
              A(j,l) = A(j,l) + as*A(k,l)
            ENDDO
            Rows(j) = Rows(j) - A(j,k)**2
          ENDDO
        ENDIF
      ELSE
        !
        !        RANK DEFICIENT PROBLEM
        Iflag = 3
        Irank = k - 1
        CALL XERMSG('SLATEC','ORTHOR',&
          'RANK OF MATRIX IS LESS THAN THE NUMBER OF ROWS.',1,1)
        RETURN
      ENDIF
    ENDDO
    GOTO 99999
  ENDIF
  Iflag = 2
  CALL XERMSG('SLATEC','ORTHOR','INVALID INPUT PARAMETERS.',2,1)
  RETURN
  !
  !
  99999 CONTINUE
  END SUBROUTINE ORTHOR
