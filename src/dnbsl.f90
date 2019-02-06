!*==DNBSL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DNBSL
      SUBROUTINE DNBSL(Abe,Lda,N,Ml,Mu,Ipvt,B,Job)
      IMPLICIT NONE
!*--DNBSL5
!***BEGIN PROLOGUE  DNBSL
!***PURPOSE  Solve a real band system using the factors computed by
!            DNBCO or DNBFA.
!***LIBRARY   SLATEC
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SNBSL-S, DNBSL-D, CNBSL-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC, SOLVE
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!     DNBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DNBCO or DNBFA.
!
!     On Entry
!
!        ABE     DOUBLE PRECISION(LDA, NC)
!                the output from DNBCO or DNBFA.
!                NC must be .GE. 2*ML+MU+1 .
!
!        LDA     INTEGER
!                the leading dimension of the array  ABE .
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
!                the pivot vector from DNBCO or DNBFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B .
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
!        setting of LDA.  It will not occur if the subroutines are
!        called correctly and if DNBCO has set RCOND .GT. 0.0
!        or DNBFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!             CALL DNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   800728  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNBSL
      INTEGER Lda , N , Ml , Mu , Ipvt(*) , Job
      DOUBLE PRECISION Abe(Lda,*) , B(*)
!
      DOUBLE PRECISION DDOT , t
      INTEGER k , kb , l , lb , ldb , lm , m , mlm , nm1
!***FIRST EXECUTABLE STATEMENT  DNBSL
      m = Mu + Ml + 1
      nm1 = N - 1
      ldb = 1 - Lda
      IF ( Job/=0 ) THEN
!
!       JOB = NONZERO, SOLVE TRANS(A) * X = B
!       FIRST SOLVE  TRANS(U)*Y = B
!
        DO k = 1 , N
          lm = MIN(k,m) - 1
          lb = k - lm
          t = DDOT(lm,Abe(k-1,Ml+2),ldb,B(lb),1)
          B(k) = (B(k)-t)/Abe(k,Ml+1)
        ENDDO
!
!       NOW SOLVE TRANS(L)*X = Y
!
        IF ( Ml/=0 ) THEN
          IF ( nm1>=1 ) THEN
            DO kb = 1 , nm1
              k = N - kb
              lm = MIN(Ml,N-k)
              mlm = Ml - (lm-1)
              B(k) = B(k) + DDOT(lm,Abe(k+lm,mlm),ldb,B(k+1),1)
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
!       JOB = 0 , SOLVE  A * X = B
!       FIRST SOLVE L*Y = B
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
              mlm = Ml - (lm-1)
              CALL DAXPY(lm,t,Abe(k+lm,mlm),ldb,B(k+1),1)
            ENDDO
          ENDIF
        ENDIF
!
!       NOW SOLVE  U*X = Y
!
        DO kb = 1 , N
          k = N + 1 - kb
          B(k) = B(k)/Abe(k,Ml+1)
          lm = MIN(k,m) - 1
          lb = k - lm
          t = -B(k)
          CALL DAXPY(lm,t,Abe(k-1,Ml+2),ldb,B(lb),1)
        ENDDO
      ENDIF
      END SUBROUTINE DNBSL
