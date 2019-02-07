!*==CGESL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CGESL
SUBROUTINE CGESL(A,Lda,N,Ipvt,B,Job)
  IMPLICIT NONE
  !*--CGESL5
  !***BEGIN PROLOGUE  CGESL
  !***PURPOSE  Solve the complex system A*X=B or CTRANS(A)*X=B using the
  !            factors computed by CGECO or CGEFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2C1
  !***TYPE      COMPLEX (SGESL-S, DGESL-D, CGESL-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     CGESL solves the complex system
  !     A * X = B  or  CTRANS(A) * X = B
  !     using the factors computed by CGECO or CGEFA.
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA, N)
  !                the output from CGECO or CGEFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from CGECO or CGEFA.
  !
  !        B       COMPLEX(N)
  !                the right hand side vector.
  !
  !        JOB     INTEGER
  !                = 0         to solve  A*X = B ,
  !                = nonzero   to solve  CTRANS(A)*X = B  where
  !                            CTRANS(A)  is the conjugate transpose.
  !
  !     On Return
  !
  !        B       the solution vector  X .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains a
  !        zero on the diagonal.  Technically this indicates singularity
  !        but it is often caused by improper arguments or improper
  !        setting of LDA .  It will not occur if the subroutines are
  !        called correctly and if CGECO has set RCOND .GT. 0.0
  !        or CGEFA has set INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL CGECO(A,LDA,N,IPVT,RCOND,Z)
  !           IF (RCOND is too small) GO TO ...
  !           DO 10 J = 1, P
  !              CALL CGESL(A,LDA,N,IPVT,C(1,J),0)
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
  !***END PROLOGUE  CGESL
  INTEGER Lda , N , Ipvt(*) , Job
  COMPLEX A(Lda,*) , B(*)
  !
  COMPLEX CDOTC , t
  INTEGER k , kb , l , nm1
  !***FIRST EXECUTABLE STATEMENT  CGESL
  nm1 = N - 1
  IF ( Job/=0 ) THEN
    !
    !        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
    !        FIRST SOLVE  CTRANS(U)*Y = B
    !
    DO k = 1 , N
      t = CDOTC(k-1,A(1,k),1,B(1),1)
      B(k) = (B(k)-t)/CONJG(A(k,k))
    ENDDO
    !
    !        NOW SOLVE CTRANS(L)*X = Y
    !
    IF ( nm1>=1 ) THEN
      DO kb = 1 , nm1
        k = N - kb
        B(k) = B(k) + CDOTC(N-k,A(k+1,k),1,B(k+1),1)
        l = Ipvt(k)
        IF ( l/=k ) THEN
          t = B(l)
          B(l) = B(k)
          B(k) = t
        ENDIF
      ENDDO
    ENDIF
  ELSE
    !
    !        JOB = 0 , SOLVE  A * X = B
    !        FIRST SOLVE  L*Y = B
    !
    IF ( nm1>=1 ) THEN
      DO k = 1 , nm1
        l = Ipvt(k)
        t = B(l)
        IF ( l/=k ) THEN
          B(l) = B(k)
          B(k) = t
        ENDIF
        CALL CAXPY(N-k,t,A(k+1,k),1,B(k+1),1)
      ENDDO
    ENDIF
    !
    !        NOW SOLVE  U*X = Y
    !
    DO kb = 1 , N
      k = N + 1 - kb
      B(k) = B(k)/A(k,k)
      t = -B(k)
      CALL CAXPY(k-1,t,A(1,k),1,B(1),1)
    ENDDO
  ENDIF
END SUBROUTINE CGESL
