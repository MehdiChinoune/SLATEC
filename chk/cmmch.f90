!*==CMMCH.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CMMCH
      SUBROUTINE CMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,Ldc,Ct,G,
     &                 Cc,Ldcc,Eps,Err,Ftl,Nout,Mv,Kprint)
      IMPLICIT NONE
!*--CMMCH6
!***BEGIN PROLOGUE  CMMCH
!***SUBSIDIARY
!***PURPOSE  Check the results of the computational tests.
!***LIBRARY   SLATEC (BLAS)
!***AUTHOR  Dongarra, J. J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!***DESCRIPTION
!
!  Checks the results of the computational tests.
!
!  Auxiliary routine for test program for Level 3 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  CMMCH
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RZERO , RONE
      PARAMETER (RZERO=0.0,RONE=1.0)
!     .. Scalar Arguments ..
      LOGICAL Ftl
      COMPLEX Alpha , Beta
      REAL Eps , Err
      INTEGER Kk , Kprint , Lda , Ldb , Ldc , Ldcc , M , N , Nout
      LOGICAL Mv
      CHARACTER*1 Transa , Transb
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldc,*) , Cc(Ldcc,*) , Ct(*)
      REAL G(*)
!     .. Local Scalars ..
      COMPLEX cl
      REAL erri
      INTEGER i , j , k
      LOGICAL ctrana , ctranb , trana , tranb
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CONJG , MAX , REAL , SQRT
!     .. Statement Functions ..
      REAL ABS1
!     .. Statement Function definitions ..
      ABS1(cl) = ABS(REAL(cl)) + ABS(AIMAG(cl))
!***FIRST EXECUTABLE STATEMENT  CMMCH
      trana = Transa=='T' .OR. Transa=='C'
      tranb = Transb=='T' .OR. Transb=='C'
      ctrana = Transa=='C'
      ctranb = Transb=='C'
!
!     Compute expected result, one column at a time, in CT using data
!     in A, B and C.
!     Compute gauges in G.
!
      DO j = 1 , N
!
        DO i = 1 , M
          Ct(i) = ZERO
          G(i) = RZERO
        ENDDO
        IF ( .NOT.trana.AND..NOT.tranb ) THEN
          DO k = 1 , Kk
            DO i = 1 , M
              Ct(i) = Ct(i) + A(i,k)*B(k,j)
              G(i) = G(i) + ABS1(A(i,k))*ABS1(B(k,j))
            ENDDO
          ENDDO
        ELSEIF ( trana.AND..NOT.tranb ) THEN
          IF ( ctrana ) THEN
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + CONJG(A(k,i))*B(k,j)
                G(i) = G(i) + ABS1(A(k,i))*ABS1(B(k,j))
              ENDDO
            ENDDO
          ELSE
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + A(k,i)*B(k,j)
                G(i) = G(i) + ABS1(A(k,i))*ABS1(B(k,j))
              ENDDO
            ENDDO
          ENDIF
        ELSEIF ( .NOT.trana.AND.tranb ) THEN
          IF ( ctranb ) THEN
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + A(i,k)*CONJG(B(j,k))
                G(i) = G(i) + ABS1(A(i,k))*ABS1(B(j,k))
              ENDDO
            ENDDO
          ELSE
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + A(i,k)*B(j,k)
                G(i) = G(i) + ABS1(A(i,k))*ABS1(B(j,k))
              ENDDO
            ENDDO
          ENDIF
        ELSEIF ( trana.AND.tranb ) THEN
          IF ( ctrana ) THEN
            IF ( ctranb ) THEN
              DO k = 1 , Kk
                DO i = 1 , M
                  Ct(i) = Ct(i) + CONJG(A(k,i))*CONJG(B(j,k))
                  G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                ENDDO
              ENDDO
            ELSE
              DO k = 1 , Kk
                DO i = 1 , M
                  Ct(i) = Ct(i) + CONJG(A(k,i))*B(j,k)
                  G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                ENDDO
              ENDDO
            ENDIF
          ELSEIF ( ctranb ) THEN
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + A(k,i)*CONJG(B(j,k))
                G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
              ENDDO
            ENDDO
          ELSE
            DO k = 1 , Kk
              DO i = 1 , M
                Ct(i) = Ct(i) + A(k,i)*B(j,k)
                G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        DO i = 1 , M
          Ct(i) = Alpha*Ct(i) + Beta*C(i,j)
          G(i) = ABS1(Alpha)*G(i) + ABS1(Beta)*ABS1(C(i,j))
        ENDDO
!
!        Compute the error ratio for this result.
!
        Err = ZERO
        DO i = 1 , M
          erri = ABS1(Ct(i)-Cc(i,j))/Eps
          IF ( G(i)/=RZERO ) erri = erri/G(i)
          Err = MAX(Err,erri)
          IF ( Err*SQRT(Eps)>=RONE ) THEN
            Ftl = .TRUE.
            IF ( Kprint>=2 ) THEN
              WRITE (Nout,FMT=99001)
              DO k = 1 , M
                IF ( Mv ) THEN
                  WRITE (Nout,FMT=99002) k , Ct(k) , Cc(k,j)
                ELSE
                  WRITE (Nout,FMT=99002) k , Cc(k,j) , Ct(k)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      RETURN
!
99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',
     &        'F ACCURATE *******',/'                       EXPECTED RE',
     &        'SULT                    COMPUTED RESULT')
99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
!
!     End of CMMCH.
!
      END SUBROUTINE CMMCH
