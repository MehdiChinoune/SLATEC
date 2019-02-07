!*==DPODI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DPODI
SUBROUTINE DPODI(A,Lda,N,Det,Job)
  IMPLICIT NONE
  !*--DPODI5
  !***BEGIN PROLOGUE  DPODI
  !***PURPOSE  Compute the determinant and inverse of a certain real
  !            symmetric positive definite matrix using the factors
  !            computed by DPOCO, DPOFA or DQRDC.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2B1B, D3B1B
  !***TYPE      DOUBLE PRECISION (SPODI-S, DPODI-D, CPODI-C)
  !***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             POSITIVE DEFINITE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DPODI computes the determinant and inverse of a certain
  !     double precision symmetric positive definite matrix (see below)
  !     using the factors computed by DPOCO, DPOFA or DQRDC.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
  !                the output  A  from DPOCO or DPOFA
  !                or the output  X  from DQRDC.
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
  !        A       If DPOCO or DPOFA was used to factor  A, then
  !                DPODI produces the upper half of INVERSE(A) .
  !                If DQRDC was used to decompose  X, then
  !                DPODI produces the upper half of inverse(TRANS(X)*X)
  !                where TRANS(X) is the transpose.
  !                Elements of  A  below the diagonal are unchanged.
  !                If the units digit of JOB is zero,  A  is unchanged.
  !
  !        DET     DOUBLE PRECISION(2)
  !                determinant of  A  or of  TRANS(X)*X  if requested.
  !                Otherwise not referenced.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. DET(1) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !     Error Condition
  !
  !        A division by zero will occur if the input factor contains
  !        a zero on the diagonal and the inverse is requested.
  !        It will not occur if the subroutines are called correctly
  !        and if DPOCO or DPOFA has set INFO .EQ. 0 .
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  DAXPY, DSCAL
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DPODI
  INTEGER Lda, N, Job
  REAL(8) :: A(Lda,*)
  REAL(8) :: Det(2)
  !
  REAL(8) :: t
  REAL(8) :: s
  INTEGER i, j, jm1, k, kp1
  !***FIRST EXECUTABLE STATEMENT  DPODI
  !
  !     COMPUTE DETERMINANT
  !
  IF ( Job/10/=0 ) THEN
    Det(1) = 1.0D0
    Det(2) = 0.0D0
    s = 10.0D0
    DO i = 1, N
      Det(1) = A(i,i)**2*Det(1)
      IF ( Det(1)==0.0D0 ) EXIT
      DO WHILE ( Det(1)<1.0D0 )
        Det(1) = s*Det(1)
        Det(2) = Det(2) - 1.0D0
      ENDDO
      DO WHILE ( Det(1)>=s )
        Det(1) = Det(1)/s
        Det(2) = Det(2) + 1.0D0
      ENDDO
    ENDDO
  ENDIF
  !
  !     COMPUTE INVERSE(R)
  !
  IF ( MOD(Job,10)/=0 ) THEN
    DO k = 1, N
      A(k,k) = 1.0D0/A(k,k)
      t = -A(k,k)
      CALL DSCAL(k-1,t,A(1,k),1)
      kp1 = k + 1
      IF ( N>=kp1 ) THEN
        DO j = kp1, N
          t = A(k,j)
          A(k,j) = 0.0D0
          CALL DAXPY(k,t,A(1,k),1,A(1,j),1)
        ENDDO
      ENDIF
    ENDDO
    !
    !        FORM  INVERSE(R) * TRANS(INVERSE(R))
    !
    DO j = 1, N
      jm1 = j - 1
      IF ( jm1>=1 ) THEN
        DO k = 1, jm1
          t = A(k,j)
          CALL DAXPY(k,t,A(1,j),1,A(1,k),1)
        ENDDO
      ENDIF
      t = A(j,j)
      CALL DSCAL(j,t,A(1,j),1)
    ENDDO
  ENDIF
END SUBROUTINE DPODI
