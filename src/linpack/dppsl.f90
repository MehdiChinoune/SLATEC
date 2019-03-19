!** DPPSL
SUBROUTINE DPPSL(Ap,N,B)
  IMPLICIT NONE
  !>
  !***
  !  Solve the real symmetric positive definite system using
  !            the factors computed by DPPCO or DPPFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B
  !***
  ! **Type:**      DOUBLE PRECISION (SPPSL-S, DPPSL-D, CPPSL-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED,
  !             POSITIVE DEFINITE, SOLVE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPPSL solves the double precision symmetric positive definite
  !     system A * X = B
  !     using the factors computed by DPPCO or DPPFA.
  !
  !     On Entry
  !
  !        AP      DOUBLE PRECISION (N*(N+1)/2)
  !                the output from DPPCO or DPPFA.
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
  !           CALL DPPCO(AP,N,RCOND,Z,INFO)
  !           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DPPSL(AP,N,C(1,J))
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
  
  INTEGER N
  REAL(8) :: Ap(*), B(*)
  !
  REAL(8) :: DDOT, t
  INTEGER k, kb, kk
  !* FIRST EXECUTABLE STATEMENT  DPPSL
  kk = 0
  DO k = 1, N
    t = DDOT(k-1,Ap(kk+1),1,B(1),1)
    kk = kk + k
    B(k) = (B(k)-t)/Ap(kk)
  ENDDO
  DO kb = 1, N
    k = N + 1 - kb
    B(k) = B(k)/Ap(kk)
    kk = kk - k
    t = -B(k)
    CALL DAXPY(k-1,t,Ap(kk+1),1,B(1),1)
  ENDDO
END SUBROUTINE DPPSL
