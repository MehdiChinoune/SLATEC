!** CHISL
SUBROUTINE CHISL(A,Lda,N,Kpvt,B)
  IMPLICIT NONE
  !>
  !***
  !  Solve the complex Hermitian system using factors obtained
  !            from CHIFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2D1A
  !***
  ! **Type:**      COMPLEX (SSISL-S, DSISL-D, CHISL-C, CSISL-C)
  !***
  ! **Keywords:**  HERMITIAN, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***
  ! **Author:**  Bunch, J., (UCSD)
  !***
  ! **Description:**
  !
  !     CHISL solves the complex Hermitian system
  !     A * X = B
  !     using the factors computed by CHIFA.
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA,N)
  !                the output from CHIFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        KVPT    INTEGER(N)
  !                the pivot vector from CHIFA.
  !
  !        B       COMPLEX(N)
  !                the right hand side vector.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero may occur if  CHICO  has set RCOND .EQ. 0.0
  !        or  CHIFA  has set INFO .NE. 0  .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL CHIFA(A,LDA,N,KVPT,INFO)
  !           IF (INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, p
  !              CALL CHISL(A,LDA,N,KVPT,C(1,J))
  !        10 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CAXPY, CDOTC

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Modified routine equivalence list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER Lda, N, Kpvt(*)
  COMPLEX A(Lda,*), B(*)
  !
  COMPLEX ak, akm1, bk, bkm1, CDOTC, denom, temp
  INTEGER k, kp
  !
  !     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
  !     D INVERSE TO B.
  !
  !* FIRST EXECUTABLE STATEMENT  CHISL
  k = N
  DO WHILE ( k/=0 )
    IF ( Kpvt(k)<0 ) THEN
      !
      !           2 X 2 PIVOT BLOCK.
      !
      IF ( k/=2 ) THEN
        kp = ABS(Kpvt(k))
        IF ( kp/=k-1 ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k-1)
          B(k-1) = B(kp)
          B(kp) = temp
        ENDIF
        !
        !              APPLY THE TRANSFORMATION.
        !
        CALL CAXPY(k-2,B(k),A(1,k),1,B(1),1)
        CALL CAXPY(k-2,B(k-1),A(1,k-1),1,B(1),1)
      ENDIF
      !
      !           APPLY D INVERSE.
      !
      ak = A(k,k)/CONJG(A(k-1,k))
      akm1 = A(k-1,k-1)/A(k-1,k)
      bk = B(k)/CONJG(A(k-1,k))
      bkm1 = B(k-1)/A(k-1,k)
      denom = ak*akm1 - 1.0E0
      B(k) = (akm1*bk-bkm1)/denom
      B(k-1) = (ak*bkm1-bk)/denom
      k = k - 2
    ELSE
      !
      !           1 X 1 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        kp = Kpvt(k)
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        ENDIF
        !
        !              APPLY THE TRANSFORMATION.
        !
        CALL CAXPY(k-1,B(k),A(1,k),1,B(1),1)
      ENDIF
      !
      !           APPLY D INVERSE.
      !
      B(k) = B(k)/A(k,k)
      k = k - 1
    ENDIF
  ENDDO
  !
  !     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
  !
  k = 1
  DO WHILE ( k<=N )
    IF ( Kpvt(k)<0 ) THEN
      !
      !           2 X 2 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        !
        !              APPLY THE TRANSFORMATION.
        !
        B(k) = B(k) + CDOTC(k-1,A(1,k),1,B(1),1)
        B(k+1) = B(k+1) + CDOTC(k-1,A(1,k+1),1,B(1),1)
        kp = ABS(Kpvt(k))
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        ENDIF
      ENDIF
      k = k + 2
    ELSE
      !
      !           1 X 1 PIVOT BLOCK.
      !
      IF ( k/=1 ) THEN
        !
        !              APPLY THE TRANSFORMATION.
        !
        B(k) = B(k) + CDOTC(k-1,A(1,k),1,B(1),1)
        kp = Kpvt(k)
        IF ( kp/=k ) THEN
          !
          !                 INTERCHANGE.
          !
          temp = B(k)
          B(k) = B(kp)
          B(kp) = temp
        ENDIF
      ENDIF
      k = k + 1
    ENDIF
  ENDDO
END SUBROUTINE CHISL
