!** DGEDI
SUBROUTINE DGEDI(A,Lda,N,Ipvt,Det,Work,Job)
  !>
  !  Compute the determinant and inverse of a matrix using the
  !            factors computed by DGECO or DGEFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D3A1, D2A1
  !***
  ! **Type:**      DOUBLE PRECISION (SGEDI-S, DGEDI-D, CGEDI-C)
  !***
  ! **Keywords:**  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DGEDI computes the determinant and inverse of a matrix
  !     using the factors computed by DGECO or DGEFA.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output from DGECO or DGEFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from DGECO or DGEFA.
  !
  !        WORK    DOUBLE PRECISION(N)
  !                work vector.  Contents destroyed.
  !
  !        JOB     INTEGER
  !                = 11   both determinant and inverse.
  !                = 01   inverse only.
  !                = 10   determinant only.
  !
  !     On Return
  !
  !        A       inverse of original matrix if requested.
  !                Otherwise unchanged.
  !
  !        DET     DOUBLE PRECISION(2)
  !                determinant of original matrix if requested.
  !                Otherwise not referenced.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains
  !        a zero on the diagonal and the inverse is requested.
  !        It will not occur if the subroutines are called correctly
  !        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
  !        INFO .EQ. 0 .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DSCAL, DSWAP

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Lda, N, Ipvt(*), Job
  REAL(DP) :: A(Lda,*), Det(2), Work(*)
  !
  REAL(DP) :: t
  REAL(DP) :: ten
  INTEGER i, j, k, kb, kp1, l, nm1
  !* FIRST EXECUTABLE STATEMENT  DGEDI
  !
  !     COMPUTE DETERMINANT
  !
  IF ( Job/10/=0 ) THEN
    Det(1) = 1.0D0
    Det(2) = 0.0D0
    ten = 10.0D0
    DO i = 1, N
      IF ( Ipvt(i)/=i ) Det(1) = -Det(1)
      Det(1) = A(i,i)*Det(1)
      IF ( Det(1)==0.0D0 ) EXIT
      DO WHILE ( ABS(Det(1))<1.0D0 )
        Det(1) = ten*Det(1)
        Det(2) = Det(2) - 1.0D0
      END DO
      DO WHILE ( ABS(Det(1))>=ten )
        Det(1) = Det(1)/ten
        Det(2) = Det(2) + 1.0D0
      END DO
    END DO
  END IF
  !
  !     COMPUTE INVERSE(U)
  !
  IF ( MOD(Job,10)/=0 ) THEN
    DO k = 1, N
      A(k,k) = 1.0D0/A(k,k)
      t = -A(k,k)
      CALL DSCAL(k-1,t,A(1,k),1)
      kp1 = k + 1
      IF ( N>=kp1 ) THEN
        DO j = kp1, N
          t = A(k,j)
          A(k,j) = 0.0D0
          CALL DAXPY(k,t,A(1,k),1,A(1,j),1)
        END DO
      END IF
    END DO
    !
    !        FORM INVERSE(U)*INVERSE(L)
    !
    nm1 = N - 1
    IF ( nm1>=1 ) THEN
      DO kb = 1, nm1
        k = N - kb
        kp1 = k + 1
        DO i = kp1, N
          Work(i) = A(i,k)
          A(i,k) = 0.0D0
        END DO
        DO j = kp1, N
          t = Work(j)
          CALL DAXPY(N,t,A(1,j),1,A(1,k),1)
        END DO
        l = Ipvt(k)
        IF ( l/=k ) CALL DSWAP(N,A(1,k),1,A(1,l),1)
      END DO
    END IF
  END IF
END SUBROUTINE DGEDI
