!** DTRDI
SUBROUTINE DTRDI(T,Ldt,N,Det,Job,Info)
  !>
  !  Compute the determinant and inverse of a triangular matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A3, D3A3
  !***
  ! **Type:**      DOUBLE PRECISION (STRDI-S, DTRDI-D, CTRDI-C)
  !***
  ! **Keywords:**  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             TRIANGULAR MATRIX
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DTRDI computes the determinant and inverse of a double precision
  !     triangular matrix.
  !
  !     On Entry
  !
  !        T       DOUBLE PRECISION(LDT,N)
  !                T contains the triangular matrix.  The zero
  !                elements of the matrix are not referenced, and
  !                the corresponding elements of the array can be
  !                used to store other information.
  !
  !        LDT     INTEGER
  !                LDT is the leading dimension of the array T.
  !
  !        N       INTEGER
  !                N is the order of the system.
  !
  !        JOB     INTEGER
  !                = 010       no det, inverse of lower triangular.
  !                = 011       no det, inverse of upper triangular.
  !                = 100       det, no inverse.
  !                = 110       det, inverse of lower triangular.
  !                = 111       det, inverse of upper triangular.
  !
  !     On Return
  !
  !        T       inverse of original matrix if requested.
  !                Otherwise unchanged.
  !
  !        DET     DOUBLE PRECISION(2)
  !                determinant of original matrix if requested.
  !                Otherwise not referenced.
  !                DETERMINANT = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !        INFO    INTEGER
  !                INFO contains zero if the system is nonsingular
  !                and the inverse is requested.
  !                Otherwise INFO contains the index of
  !                a zero diagonal element of T.
  !
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Ldt, N, Job, Info
  REAL(DP) :: T(Ldt,*), Det(2)
  !
  REAL(DP) :: temp
  REAL(DP) :: ten
  INTEGER i, j, k, kb, km1, kp1
  !* FIRST EXECUTABLE STATEMENT  DTRDI
  !
  !        COMPUTE DETERMINANT
  !
  IF ( Job/100/=0 ) THEN
    Det(1) = 1.0D0
    Det(2) = 0.0D0
    ten = 10.0D0
    DO i = 1, N
      Det(1) = T(i,i)*Det(1)
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
  !        COMPUTE INVERSE OF UPPER TRIANGULAR
  !
  IF ( MOD(Job/10,10)/=0 ) THEN
    IF ( MOD(Job,10)==0 ) THEN
      !
      !              COMPUTE INVERSE OF LOWER TRIANGULAR
      !
      DO kb = 1, N
        k = N + 1 - kb
        Info = k
        IF ( T(k,k)==0.0D0 ) RETURN
        T(k,k) = 1.0D0/T(k,k)
        temp = -T(k,k)
        IF ( k/=N ) CALL DSCAL(N-k,temp,T(k+1,k),1)
        km1 = k - 1
        IF ( km1>=1 ) THEN
          DO j = 1, km1
            temp = T(k,j)
            T(k,j) = 0.0D0
            CALL DAXPY(N-k+1,temp,T(k,k),1,T(k,j),1)
          END DO
        END IF
      END DO
      Info = 0
    ELSE
      DO k = 1, N
        Info = k
        IF ( T(k,k)==0.0D0 ) RETURN
        T(k,k) = 1.0D0/T(k,k)
        temp = -T(k,k)
        CALL DSCAL(k-1,temp,T(1,k),1)
        kp1 = k + 1
        IF ( N>=kp1 ) THEN
          DO j = kp1, N
            temp = T(k,j)
            T(k,j) = 0.0D0
            CALL DAXPY(k,temp,T(1,k),1,T(1,j),1)
          END DO
        END IF
      END DO
      Info = 0
    END IF
  END IF
  RETURN
END SUBROUTINE DTRDI
