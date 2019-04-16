!** DPOSL
SUBROUTINE DPOSL(A,Lda,N,B)
  !>
  !***
  !  Solve the real symmetric positive definite linear system
  !            using the factors computed by DPOCO or DPOFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B
  !***
  ! **Type:**      DOUBLE PRECISION (SPOSL-S, DPOSL-D, CPOSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPOSL solves the double precision symmetric positive definite
  !     system A * X = B
  !     using the factors computed by DPOCO or DPOFA.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output from DPOCO or DPOFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        B       DOUBLE PRECISION(N)
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
  !        singularity, but it is usually caused by improper subroutine
  !        arguments.  It will not occur if the subroutines are called
  !        correctly and  INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
  !           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DPOSL(A,LDA,N,C(1,J))
  !        10 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Lda, N
  REAL(8) :: A(Lda,*), B(*)
  !
  REAL(8) :: t
  INTEGER k, kb
  !
  !     SOLVE TRANS(R)*Y = B
  !
  !* FIRST EXECUTABLE STATEMENT  DPOSL
  DO k = 1, N
    t = DDOT(k-1,A(1,k),1,B(1),1)
    B(k) = (B(k)-t)/A(k,k)
  END DO
  !
  !     SOLVE R*X = Y
  !
  DO kb = 1, N
    k = N + 1 - kb
    B(k) = B(k)/A(k,k)
    t = -B(k)
    CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
  END DO
END SUBROUTINE DPOSL
