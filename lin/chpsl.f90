!*==CHPSL.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CHPSL
      SUBROUTINE CHPSL(Ap,N,Kpvt,B)
      IMPLICIT NONE
!*--CHPSL5
!***BEGIN PROLOGUE  CHPSL
!***PURPOSE  Solve a complex Hermitian system using factors obtained
!            from CHPFA.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1A
!***TYPE      COMPLEX (SSPSL-S, DSPSL-D, CHPSL-C, CSPSL-C)
!***KEYWORDS  HERMITIAN, LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, SOLVE
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     CHISL solves the complex Hermitian system
!     A * X = B
!     using the factors computed by CHPFA.
!
!     On Entry
!
!        AP      COMPLEX(N*(N+1)/2)
!                the output from CHPFA.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        KVPT    INTEGER(N)
!                the pivot vector from CHPFA.
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
!        A division by zero may occur if  CHPCO  has set RCOND .EQ. 0.0
!        or  CHPFA  has set INFO .NE. 0  .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL CHPFA(AP,N,KVPT,INFO)
!           IF (INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL CHPSL(AP,N,KVPT,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Modified routine equivalence list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CHPSL
      INTEGER N , Kpvt(*)
      COMPLEX Ap(*) , B(*)
!
      COMPLEX ak , akm1 , bk , bkm1 , CDOTC , denom , temp
      INTEGER ik , ikm1 , ikp1 , k , kk , km1k , km1km1 , kp
!
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
!     D INVERSE TO B.
!
!***FIRST EXECUTABLE STATEMENT  CHPSL
      k = N
      ik = (N*(N-1))/2
      DO WHILE ( k/=0 )
        kk = ik + k
        IF ( Kpvt(k)<0 ) THEN
!
!           2 X 2 PIVOT BLOCK.
!
          ikm1 = ik - (k-1)
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
            CALL CAXPY(k-2,B(k),Ap(ik+1),1,B(1),1)
            CALL CAXPY(k-2,B(k-1),Ap(ikm1+1),1,B(1),1)
          ENDIF
!
!           APPLY D INVERSE.
!
          km1k = ik + k - 1
          kk = ik + k
          ak = Ap(kk)/CONJG(Ap(km1k))
          km1km1 = ikm1 + k - 1
          akm1 = Ap(km1km1)/Ap(km1k)
          bk = B(k)/CONJG(Ap(km1k))
          bkm1 = B(k-1)/Ap(km1k)
          denom = ak*akm1 - 1.0E0
          B(k) = (akm1*bk-bkm1)/denom
          B(k-1) = (ak*bkm1-bk)/denom
          k = k - 2
          ik = ik - (k+1) - k
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
            CALL CAXPY(k-1,B(k),Ap(ik+1),1,B(1),1)
          ENDIF
!
!           APPLY D INVERSE.
!
          B(k) = B(k)/Ap(kk)
          k = k - 1
          ik = ik - k
        ENDIF
      ENDDO
!
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
!
      k = 1
      ik = 0
      DO WHILE ( k<=N )
        IF ( Kpvt(k)<0 ) THEN
!
!           2 X 2 PIVOT BLOCK.
!
          IF ( k/=1 ) THEN
!
!              APPLY THE TRANSFORMATION.
!
            B(k) = B(k) + CDOTC(k-1,Ap(ik+1),1,B(1),1)
            ikp1 = ik + k
            B(k+1) = B(k+1) + CDOTC(k-1,Ap(ikp1+1),1,B(1),1)
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
          ik = ik + k + k + 1
          k = k + 2
        ELSE
!
!           1 X 1 PIVOT BLOCK.
!
          IF ( k/=1 ) THEN
!
!              APPLY THE TRANSFORMATION.
!
            B(k) = B(k) + CDOTC(k-1,Ap(ik+1),1,B(1),1)
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
          ik = ik + k
          k = k + 1
        ENDIF
      ENDDO
      END SUBROUTINE CHPSL
