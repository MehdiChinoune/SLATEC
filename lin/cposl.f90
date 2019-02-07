!*==CPOSL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CPOSL
SUBROUTINE CPOSL(A,Lda,N,B)
  IMPLICIT NONE
  !*--CPOSL5
  !***BEGIN PROLOGUE  CPOSL
  !***PURPOSE  Solve the complex Hermitian positive definite linear system
  !            using the factors computed by CPOCO or CPOFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D1B
  !***TYPE      COMPLEX (SPOSL-S, DPOSL-D, CPOSL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     CPOSL solves the COMPLEX Hermitian positive definite system
  !     A * X = B
  !     using the factors computed by CPOCO or CPOFA.
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA, N)
  !                the output from CPOCO or CPOFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
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
  !           CALL CPOCO(A,LDA,N,RCOND,Z,INFO)
  !           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL CPOSL(A,LDA,N,C(1,J))
  !        10 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CAXPY, CDOTC
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CPOSL
  INTEGER Lda , N
  COMPLEX A(Lda,*) , B(*)
  !
  COMPLEX CDOTC , t
  INTEGER k , kb
  !
  !     SOLVE CTRANS(R)*Y = B
  !
  !***FIRST EXECUTABLE STATEMENT  CPOSL
  DO k = 1 , N
    t = CDOTC(k-1,A(1,k),1,B(1),1)
    B(k) = (B(k)-t)/A(k,k)
  ENDDO
  !
  !     SOLVE R*X = Y
  !
  DO kb = 1 , N
    k = N + 1 - kb
    B(k) = B(k)/A(k,k)
    t = -B(k)
    CALL CAXPY(k-1,t,A(1,k),1,B(1),1)
  ENDDO
END SUBROUTINE CPOSL
