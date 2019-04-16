!** CPBSL
SUBROUTINE CPBSL(Abd,Lda,N,M,B)
  !>
  !***
  !  Solve the complex Hermitian positive definite band system
  !            using the factors computed by CPBCO or CPBFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2D2
  !***
  ! **Type:**      COMPLEX (SPBSL-S, DPBSL-D, CPBSL-C)
  !***
  ! **Keywords:**  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             POSITIVE DEFINITE, SOLVE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     CPBSL solves the complex Hermitian positive definite band
  !     system  A*X = B
  !     using the factors computed by CPBCO or CPBFA.
  !
  !     On Entry
  !
  !        ABD     COMPLEX(LDA, N)
  !                the output from CPBCO or CPBFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        M       INTEGER
  !                the number of diagonals above the main diagonal.
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
  !        A division by zero will occur if the input factor contains
  !        a zero on the diagonal.  Technically this indicates
  !        singularity but it is usually caused by improper subroutine
  !        arguments.  It will not occur if the subroutines are called
  !        correctly and  INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL CPBCO(ABD,LDA,N,RCOND,Z,INFO)
  !           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL CPBSL(ABD,LDA,N,C(1,J))
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
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Lda, N, M
  COMPLEX Abd(Lda,*), B(*)
  !
  COMPLEX t
  INTEGER k, kb, la, lb, lm
  !
  !     SOLVE CTRANS(R)*Y = B
  !
  !* FIRST EXECUTABLE STATEMENT  CPBSL
  DO k = 1, N
    lm = MIN(k-1,M)
    la = M + 1 - lm
    lb = k - lm
    t = CDOTC(lm,Abd(la,k),1,B(lb),1)
    B(k) = (B(k)-t)/Abd(M+1,k)
  END DO
  !
  !     SOLVE R*X = Y
  !
  DO kb = 1, N
    k = N + 1 - kb
    lm = MIN(k-1,M)
    la = M + 1 - lm
    lb = k - lm
    B(k) = B(k)/Abd(M+1,k)
    t = -B(k)
    CALL CAXPY(lm,t,Abd(la,k),1,B(lb),1)
  END DO
END SUBROUTINE CPBSL
