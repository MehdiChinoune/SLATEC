!*==SMMCH.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK SMMCH
      SUBROUTINE SMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,Ldc,Ct,G,
     &                 Cc,Ldcc,Eps,Err,Ftl,Nout,Mv,Kprint)
      IMPLICIT NONE
!*--SMMCH6
!***BEGIN PROLOGUE  SMMCH
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
!***END PROLOGUE  SMMCH
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
!     .. Scalar Arguments ..
      REAL Alpha , Beta , Eps , Err
      INTEGER Kk , Kprint , Lda , Ldb , Ldc , Ldcc , M , N , Nout
      LOGICAL Mv , Ftl
      CHARACTER*1 Transa , Transb
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*) , C(Ldc,*) , Cc(Ldcc,*) , Ct(*) , G(*)
!     .. Local Scalars ..
      REAL erri
      INTEGER i , j , k
      LOGICAL trana , tranb
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!***FIRST EXECUTABLE STATEMENT  SMMCH
      trana = Transa=='T' .OR. Transa=='C'
      tranb = Transb=='T' .OR. Transb=='C'
!
!     Compute expected result, one column at a time, in CT using data
!     in A, B and C.
!     Compute gauges in G.
!
      DO j = 1 , N
!
        DO i = 1 , M
          Ct(i) = ZERO
          G(i) = ZERO
        ENDDO
        IF ( .NOT.trana.AND..NOT.tranb ) THEN
          DO k = 1 , Kk
            DO i = 1 , M
              Ct(i) = Ct(i) + A(i,k)*B(k,j)
              G(i) = G(i) + ABS(A(i,k))*ABS(B(k,j))
            ENDDO
          ENDDO
        ELSEIF ( trana.AND..NOT.tranb ) THEN
          DO k = 1 , Kk
            DO i = 1 , M
              Ct(i) = Ct(i) + A(k,i)*B(k,j)
              G(i) = G(i) + ABS(A(k,i))*ABS(B(k,j))
            ENDDO
          ENDDO
        ELSEIF ( .NOT.trana.AND.tranb ) THEN
          DO k = 1 , Kk
            DO i = 1 , M
              Ct(i) = Ct(i) + A(i,k)*B(j,k)
              G(i) = G(i) + ABS(A(i,k))*ABS(B(j,k))
            ENDDO
          ENDDO
        ELSEIF ( trana.AND.tranb ) THEN
          DO k = 1 , Kk
            DO i = 1 , M
              Ct(i) = Ct(i) + A(k,i)*B(j,k)
              G(i) = G(i) + ABS(A(k,i))*ABS(B(j,k))
            ENDDO
          ENDDO
        ENDIF
        DO i = 1 , M
          Ct(i) = Alpha*Ct(i) + Beta*C(i,j)
          G(i) = ABS(Alpha)*G(i) + ABS(Beta)*ABS(C(i,j))
        ENDDO
!
!        Compute the error ratio for this result.
!
        Err = ZERO
        DO i = 1 , M
          erri = ABS(Ct(i)-Cc(i,j))/Eps
          IF ( G(i)/=ZERO ) erri = erri/G(i)
          Err = MAX(Err,erri)
          IF ( Err*SQRT(Eps)>=ONE ) THEN
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
     &        'F ACCURATE *******',/'           EXPECTED RESULT   COMPU',
     &        'TED RESULT')
99002 FORMAT (1X,I7,2G18.6)
!
!     End of SMMCH.
!
      END SUBROUTINE SMMCH
