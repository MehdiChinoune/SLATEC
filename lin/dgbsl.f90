!*==DGBSL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DGBSL
SUBROUTINE DGBSL(Abd,Lda,N,Ml,Mu,Ipvt,B,Job)
  IMPLICIT NONE
  !*--DGBSL5
  !***BEGIN PROLOGUE  DGBSL
  !***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
  !            the factors computed by DGBCO or DGBFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2A2
  !***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
  !***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGBSL solves the double precision band system
  !     A * X = B  or  TRANS(A) * X = B
  !     using the factors computed by DGBCO or DGBFA.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the output from DGBCO or DGBFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from DGBCO or DGBFA.
  !
  !        B       DOUBLE PRECISION(N)
  !                the right hand side vector.
  !
  !        JOB     INTEGER
  !                = 0         to solve  A*X = B ,
  !                = nonzero   to solve  TRANS(A)*X = B , where
  !                            TRANS(A)  is the transpose.
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
  !        called correctly and if DGBCO has set RCOND .GT. 0.0
  !        or DGBFA has set INFO .EQ. 0 .
  !
  !     To compute  INVERSE(A) * C  where  C  is a matrix
  !     with  P  columns
  !           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
  !           IF (RCOND is too small) GO TO ...
  !           DO 10 J = 1, P
  !              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
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
  !***END PROLOGUE  DGBSL
  INTEGER Lda , N , Ml , Mu , Ipvt(*) , Job
  DOUBLE PRECISION Abd(Lda,*) , B(*)
  !
  DOUBLE PRECISION DDOT , t
  INTEGER k , kb , l , la , lb , lm , m , nm1
  !***FIRST EXECUTABLE STATEMENT  DGBSL
  m = Mu + Ml + 1
  nm1 = N - 1
  IF ( Job/=0 ) THEN
    !
    !        JOB = NONZERO, SOLVE  TRANS(A) * X = B
    !        FIRST SOLVE  TRANS(U)*Y = B
    !
    DO k = 1 , N
      lm = MIN(k,m) - 1
      la = m - lm
      lb = k - lm
      t = DDOT(lm,Abd(la,k),1,B(lb),1)
      B(k) = (B(k)-t)/Abd(m,k)
    ENDDO
    !
    !        NOW SOLVE TRANS(L)*X = Y
    !
    IF ( Ml/=0 ) THEN
      IF ( nm1>=1 ) THEN
        DO kb = 1 , nm1
          k = N - kb
          lm = MIN(Ml,N-k)
          B(k) = B(k) + DDOT(lm,Abd(m+1,k),1,B(k+1),1)
          l = Ipvt(k)
          IF ( l/=k ) THEN
            t = B(l)
            B(l) = B(k)
            B(k) = t
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ELSE
    !
    !        JOB = 0 , SOLVE  A * X = B
    !        FIRST SOLVE L*Y = B
    !
    IF ( Ml/=0 ) THEN
      IF ( nm1>=1 ) THEN
        DO k = 1 , nm1
          lm = MIN(Ml,N-k)
          l = Ipvt(k)
          t = B(l)
          IF ( l/=k ) THEN
            B(l) = B(k)
            B(k) = t
          ENDIF
          CALL DAXPY(lm,t,Abd(m+1,k),1,B(k+1),1)
        ENDDO
      ENDIF
    ENDIF
    !
    !        NOW SOLVE  U*X = Y
    !
    DO kb = 1 , N
      k = N + 1 - kb
      B(k) = B(k)/Abd(m,k)
      lm = MIN(k,m) - 1
      la = m - lm
      lb = k - lm
      t = -B(k)
      CALL DAXPY(lm,t,Abd(la,k),1,B(lb),1)
    ENDDO
  ENDIF
END SUBROUTINE DGBSL
