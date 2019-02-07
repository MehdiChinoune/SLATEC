!*==DPBSL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DPBSL
SUBROUTINE DPBSL(Abd,Lda,N,M,B)
  IMPLICIT NONE
  !*--DPBSL5
  !***BEGIN PROLOGUE  DPBSL
  !***PURPOSE  Solve a real symmetric positive definite band system
  !            using the factors computed by DPBCO or DPBFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2B2
  !***TYPE      DOUBLE PRECISION (SPBSL-S, DPBSL-D, CPBSL-C)
  !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             POSITIVE DEFINITE, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DPBSL solves the double precision symmetric positive definite
  !     band system  A*X = B
  !     using the factors computed by DPBCO or DPBFA.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the output from DPBCO or DPBFA.
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
  !        correctly, and  INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)
  !           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DPBSL(ABD,LDA,N,C(1,J))
  !        10 CONTINUE
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DPBSL
  INTEGER Lda , N , M
  REAL(8) :: Abd(Lda,*) , B(*)
  !
  REAL(8) :: DDOT , t
  INTEGER k , kb , la , lb , lm
  !
  !     SOLVE TRANS(R)*Y = B
  !
  !***FIRST EXECUTABLE STATEMENT  DPBSL
  DO k = 1 , N
    lm = MIN(k-1,M)
    la = M + 1 - lm
    lb = k - lm
    t = DDOT(lm,Abd(la,k),1,B(lb),1)
    B(k) = (B(k)-t)/Abd(M+1,k)
  ENDDO
  !
  !     SOLVE R*X = Y
  !
  DO kb = 1 , N
    k = N + 1 - kb
    lm = MIN(k-1,M)
    la = M + 1 - lm
    lb = k - lm
    B(k) = B(k)/Abd(M+1,k)
    t = -B(k)
    CALL DAXPY(lm,t,Abd(la,k),1,B(lb),1)
  ENDDO
END SUBROUTINE DPBSL
