!** DPPDI
SUBROUTINE DPPDI(Ap,N,Det,Job)
  !>
  !***
  !  Compute the determinant and inverse of a real symmetric
  !            positive definite matrix using factors from DPPCO or DPPFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B, D3B1B
  !***
  ! **Type:**      DOUBLE PRECISION (SPPDI-S, DPPDI-D, CPPDI-C)
  !***
  ! **Keywords:**  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             PACKED, POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPPDI computes the determinant and inverse
  !     of a double precision symmetric positive definite matrix
  !     using the factors computed by DPPCO or DPPFA .
  !
  !     On Entry
  !
  !        AP      DOUBLE PRECISION (N*(N+1)/2)
  !                the output from DPPCO or DPPFA.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        JOB     INTEGER
  !                = 11   both determinant and inverse.
  !                = 01   inverse only.
  !                = 10   determinant only.
  !
  !     On Return
  !
  !        AP      the upper triangular half of the inverse .
  !                The strict lower triangle is unaltered.
  !
  !        DET     DOUBLE PRECISION(2)
  !                determinant of original matrix if requested.
  !                Otherwise not referenced.
  !                DETERMINANT = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. DET(1) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains
  !        a zero on the diagonal and the inverse is requested.
  !        It will not occur if the subroutines are called correctly
  !        and if DPOCO or DPOFA has set INFO .EQ. 0 .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER N, Job
  REAL(8) :: Ap(*)
  REAL(8) :: Det(2)
  !
  REAL(8) :: t
  REAL(8) :: s
  INTEGER i, ii, j, jj, jm1, j1, k, kj, kk, kp1, k1
  !* FIRST EXECUTABLE STATEMENT  DPPDI
  !
  !     COMPUTE DETERMINANT
  !
  IF ( Job/10/=0 ) THEN
    Det(1) = 1.0D0
    Det(2) = 0.0D0
    s = 10.0D0
    ii = 0
    DO i = 1, N
      ii = ii + i
      Det(1) = Ap(ii)**2*Det(1)
      IF ( Det(1)==0.0D0 ) EXIT
      DO WHILE ( Det(1)<1.0D0 )
        Det(1) = s*Det(1)
        Det(2) = Det(2) - 1.0D0
      END DO
      DO WHILE ( Det(1)>=s )
        Det(1) = Det(1)/s
        Det(2) = Det(2) + 1.0D0
      END DO
    END DO
  END IF
  !
  !     COMPUTE INVERSE(R)
  !
  IF ( MOD(Job,10)/=0 ) THEN
    kk = 0
    DO k = 1, N
      k1 = kk + 1
      kk = kk + k
      Ap(kk) = 1.0D0/Ap(kk)
      t = -Ap(kk)
      CALL DSCAL(k-1,t,Ap(k1),1)
      kp1 = k + 1
      j1 = kk + 1
      kj = kk + k
      IF ( N>=kp1 ) THEN
        DO j = kp1, N
          t = Ap(kj)
          Ap(kj) = 0.0D0
          CALL DAXPY(k,t,Ap(k1),1,Ap(j1),1)
          j1 = j1 + j
          kj = kj + j
        END DO
      END IF
    END DO
    !
    !        FORM  INVERSE(R) * TRANS(INVERSE(R))
    !
    jj = 0
    DO j = 1, N
      j1 = jj + 1
      jj = jj + j
      jm1 = j - 1
      k1 = 1
      kj = j1
      IF ( jm1>=1 ) THEN
        DO k = 1, jm1
          t = Ap(kj)
          CALL DAXPY(k,t,Ap(j1),1,Ap(k1),1)
          k1 = k1 + k
          kj = kj + 1
        END DO
      END IF
      t = Ap(jj)
      CALL DSCAL(j,t,Ap(j1),1)
    END DO
  END IF
END SUBROUTINE DPPDI
