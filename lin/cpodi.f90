!*==CPODI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CPODI
      SUBROUTINE CPODI(A,Lda,N,Det,Job)
      IMPLICIT NONE
!*--CPODI5
!***BEGIN PROLOGUE  CPODI
!***PURPOSE  Compute the determinant and inverse of a certain complex
!            Hermitian positive definite matrix using the factors
!            computed by CPOCO, CPOFA, or CQRDC.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B, D3D1B
!***TYPE      COMPLEX (SPODI-S, DPODI-D, CPODI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPODI computes the determinant and inverse of a certain
!     complex Hermitian positive definite matrix (see below)
!     using the factors computed by CPOCO, CPOFA or CQRDC.
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the output  A  from CPOCO or CPOFA
!                or the output  X  from CQRDC.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        A       If CPOCO or CPOFA was used to factor  A  then
!                CPODI produces the upper half of INVERSE(A) .
!                If CQRDC was used to decompose  X  then
!                CPODI produces the upper half of INVERSE(CTRANS(X)*X)
!                where CTRANS(X) is the conjugate transpose.
!                Elements of  A  below the diagonal are unchanged.
!                If the units digit of JOB is zero,  A  is unchanged.
!
!        DET     REAL(2)
!                determinant of  A  or of  CTRANS(X)*X  if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0 .LE. DET(1) .LT. 10.0
!                or  DET(1) .EQ. 0.0 .
!
!     Error Condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if CPOCO or CPOFA has set INFO .EQ. 0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPODI
      INTEGER Lda , N , Job
      COMPLEX A(Lda,*)
      REAL Det(2)
!
      COMPLEX t
      REAL s
      INTEGER i , j , jm1 , k , kp1
!***FIRST EXECUTABLE STATEMENT  CPODI
!
!     COMPUTE DETERMINANT
!
      IF ( Job/10/=0 ) THEN
        Det(1) = 1.0E0
        Det(2) = 0.0E0
        s = 10.0E0
        DO i = 1 , N
          Det(1) = REAL(A(i,i))**2*Det(1)
          IF ( Det(1)==0.0E0 ) EXIT
          DO WHILE ( Det(1)<1.0E0 )
            Det(1) = s*Det(1)
            Det(2) = Det(2) - 1.0E0
          ENDDO
          DO WHILE ( Det(1)>=s )
            Det(1) = Det(1)/s
            Det(2) = Det(2) + 1.0E0
          ENDDO
        ENDDO
      ENDIF
!
!     COMPUTE INVERSE(R)
!
      IF ( MOD(Job,10)/=0 ) THEN
        DO k = 1 , N
          A(k,k) = (1.0E0,0.0E0)/A(k,k)
          t = -A(k,k)
          CALL CSCAL(k-1,t,A(1,k),1)
          kp1 = k + 1
          IF ( N>=kp1 ) THEN
            DO j = kp1 , N
              t = A(k,j)
              A(k,j) = (0.0E0,0.0E0)
              CALL CAXPY(k,t,A(1,k),1,A(1,j),1)
            ENDDO
          ENDIF
        ENDDO
!
!        FORM  INVERSE(R) * CTRANS(INVERSE(R))
!
        DO j = 1 , N
          jm1 = j - 1
          IF ( jm1>=1 ) THEN
            DO k = 1 , jm1
              t = CONJG(A(k,j))
              CALL CAXPY(k,t,A(1,j),1,A(1,k),1)
            ENDDO
          ENDIF
          t = CONJG(A(j,j))
          CALL CSCAL(j,t,A(1,j),1)
        ENDDO
      ENDIF
      END SUBROUTINE CPODI
