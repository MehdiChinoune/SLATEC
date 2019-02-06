!*==CSIDI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CSIDI
      SUBROUTINE CSIDI(A,Lda,N,Kpvt,Det,Work,Job)
      IMPLICIT NONE
!*--CSIDI5
!***BEGIN PROLOGUE  CSIDI
!***PURPOSE  Compute the determinant and inverse of a complex symmetric
!            matrix using the factors from CSIFA.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C1, D3C1
!***TYPE      COMPLEX (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     CSIDI computes the determinant and inverse
!     of a complex symmetric matrix using the factors from CSIFA.
!
!     On Entry
!
!        A       COMPLEX(LDA,N)
!                the output from CSIFA.
!
!        LDA     INTEGER
!                the leading dimension of the array A .
!
!        N       INTEGER
!                the order of the matrix A .
!
!        KVPT    INTEGER(N)
!                the pivot vector from CSIFA.
!
!        WORK    COMPLEX(N)
!                work vector.  Contents destroyed.
!
!        JOB     INTEGER
!                JOB has the decimal expansion  AB  where
!                   If  B .NE. 0, the inverse is computed,
!                   If  A .NE. 0, the determinant is computed,
!
!                For example, JOB = 11  gives both.
!
!     On Return
!
!        Variables not requested by JOB are not used.
!
!        A      contains the upper triangle of the inverse of
!               the original matrix.  The strict lower triangle
!               is never referenced.
!
!        DET    COMPLEX(2)
!               determinant of original matrix.
!               Determinant = DET(1) * 10.0**DET(2)
!               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
!               or DET(1) = 0.0.
!
!     Error Condition
!
!        A division by zero may occur if the inverse is requested
!        and  CSICO  has set RCOND .EQ. 0.0
!        or  CSIFA  has set  INFO .NE. 0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CCOPY, CDOTU, CSWAP
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Corrected category and modified routine equivalence
!           list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CSIDI
      INTEGER Lda , N , Job
      COMPLEX A(Lda,*) , Det(2) , Work(*)
      INTEGER Kpvt(*)
!
      COMPLEX ak , akp1 , akkp1 , CDOTU , d , t , temp
      REAL ten
      INTEGER j , jb , k , km1 , ks , kstep
      LOGICAL noinv , nodet
      REAL, EXTERNAL :: CABS1
!
!***FIRST EXECUTABLE STATEMENT  CSIDI
      noinv = MOD(Job,10)==0
      nodet = MOD(Job,100)/10==0
!
      IF ( .NOT.(nodet) ) THEN
        Det(1) = (1.0E0,0.0E0)
        Det(2) = (0.0E0,0.0E0)
        ten = 10.0E0
        t = (0.0E0,0.0E0)
        DO k = 1 , N
          d = A(k,k)
!
!           CHECK IF 1 BY 1
!
          IF ( Kpvt(k)<=0 ) THEN
!
!              2 BY 2 BLOCK
!              USE DET (D  T)  =  (D/T * C - T) * T
!                      (T  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
            IF ( CABS1(t)/=0.0E0 ) THEN
              d = t
              t = (0.0E0,0.0E0)
            ELSE
              t = A(k,k+1)
              d = (d/t)*A(k+1,k+1) - t
            ENDIF
          ENDIF
!
          Det(1) = d*Det(1)
          IF ( CABS1(Det(1))/=0.0E0 ) THEN
            DO WHILE ( CABS1(Det(1))<1.0E0 )
              Det(1) = CMPLX(ten,0.0E0)*Det(1)
              Det(2) = Det(2) - (1.0E0,0.0E0)
            ENDDO
            DO WHILE ( CABS1(Det(1))>=ten )
              Det(1) = Det(1)/CMPLX(ten,0.0E0)
              Det(2) = Det(2) + (1.0E0,0.0E0)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
!
!     COMPUTE INVERSE(A)
!
      IF ( .NOT.(noinv) ) THEN
        k = 1
        DO WHILE ( k<=N )
          km1 = k - 1
          IF ( Kpvt(k)<0 ) THEN
!
!              2 BY 2
!
            t = A(k,k+1)
            ak = A(k,k)/t
            akp1 = A(k+1,k+1)/t
            akkp1 = A(k,k+1)/t
            d = t*(ak*akp1-(1.0E0,0.0E0))
            A(k,k) = akp1/d
            A(k+1,k+1) = ak/d
            A(k,k+1) = -akkp1/d
            IF ( km1>=1 ) THEN
              CALL CCOPY(km1,A(1,k+1),1,Work,1)
              DO j = 1 , km1
                A(j,k+1) = CDOTU(j,A(1,j),1,Work,1)
                CALL CAXPY(j-1,Work(j),A(1,j),1,A(1,k+1),1)
              ENDDO
              A(k+1,k+1) = A(k+1,k+1) + CDOTU(km1,Work,1,A(1,k+1),1)
              A(k,k+1) = A(k,k+1) + CDOTU(km1,A(1,k),1,A(1,k+1),1)
              CALL CCOPY(km1,A(1,k),1,Work,1)
              DO j = 1 , km1
                A(j,k) = CDOTU(j,A(1,j),1,Work,1)
                CALL CAXPY(j-1,Work(j),A(1,j),1,A(1,k),1)
              ENDDO
              A(k,k) = A(k,k) + CDOTU(km1,Work,1,A(1,k),1)
            ENDIF
            kstep = 2
          ELSE
!
!              1 BY 1
!
            A(k,k) = (1.0E0,0.0E0)/A(k,k)
            IF ( km1>=1 ) THEN
              CALL CCOPY(km1,A(1,k),1,Work,1)
              DO j = 1 , km1
                A(j,k) = CDOTU(j,A(1,j),1,Work,1)
                CALL CAXPY(j-1,Work(j),A(1,j),1,A(1,k),1)
              ENDDO
              A(k,k) = A(k,k) + CDOTU(km1,Work,1,A(1,k),1)
            ENDIF
            kstep = 1
          ENDIF
!
!           SWAP
!
          ks = ABS(Kpvt(k))
          IF ( ks/=k ) THEN
            CALL CSWAP(ks,A(1,ks),1,A(1,k),1)
            DO jb = ks , k
              j = k + ks - jb
              temp = A(j,k)
              A(j,k) = A(ks,j)
              A(ks,j) = temp
            ENDDO
            IF ( kstep/=1 ) THEN
              temp = A(ks,k+1)
              A(ks,k+1) = A(k,k+1)
              A(k,k+1) = temp
            ENDIF
          ENDIF
          k = k + kstep
        ENDDO
      ENDIF
      END SUBROUTINE CSIDI
