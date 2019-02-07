!*==STRDI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK STRDI
SUBROUTINE STRDI(T,Ldt,N,Det,Job,Info)
  IMPLICIT NONE
  !*--STRDI5
  !***BEGIN PROLOGUE  STRDI
  !***PURPOSE  Compute the determinant and inverse of a triangular matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2A3, D3A3
  !***TYPE      SINGLE PRECISION (STRDI-S, DTRDI-D, CTRDI-C)
  !***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             TRIANGULAR
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     STRDI computes the determinant and inverse of a real
  !     triangular matrix.
  !
  !     On Entry
  !
  !        T       REAL(LDT,N)
  !                T contains the triangular matrix.  The zero
  !                elements of the matrix are not referenced, and
  !                the corresponding elements of the array can be
  !                used to store other information.
  !
  !        LDT     INTEGER
  !                LDT is the leading dimension of the array T.
  !
  !        N       INTEGER
  !                N is the order of the system.
  !
  !        JOB     INTEGER
  !                = 010       no det, inverse of lower triangular.
  !                = 011       no det, inverse of upper triangular.
  !                = 100       det, no inverse.
  !                = 110       det, inverse of lower triangular.
  !                = 111       det, inverse of upper triangular.
  !
  !     On Return
  !
  !        T       inverse of original matrix if requested.
  !                Otherwise unchanged.
  !
  !        DET     REAL(2)
  !                determinant of original matrix if requested.
  !                Otherwise not referenced.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !        INFO    INTEGER
  !                INFO contains zero if the system is nonsingular
  !                and the inverse is requested.
  !                Otherwise INFO contains the index of
  !                a zero diagonal element of T.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  SAXPY, SSCAL
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  STRDI
  INTEGER Ldt , N , Job , Info
  REAL T(Ldt,*) , Det(2)
  !
  REAL temp
  REAL ten
  INTEGER i , j , k , kb , km1 , kp1
  !***FIRST EXECUTABLE STATEMENT  STRDI
  !
  !        COMPUTE DETERMINANT
  !
  IF ( Job/100/=0 ) THEN
    Det(1) = 1.0E0
    Det(2) = 0.0E0
    ten = 10.0E0
    DO i = 1 , N
      Det(1) = T(i,i)*Det(1)
      IF ( Det(1)==0.0E0 ) EXIT
      DO WHILE ( ABS(Det(1))<1.0E0 )
        Det(1) = ten*Det(1)
        Det(2) = Det(2) - 1.0E0
      ENDDO
      DO WHILE ( ABS(Det(1))>=ten )
        Det(1) = Det(1)/ten
        Det(2) = Det(2) + 1.0E0
      ENDDO
    ENDDO
  ENDIF
  !
  !        COMPUTE INVERSE OF UPPER TRIANGULAR
  !
  IF ( MOD(Job/10,10)/=0 ) THEN
    IF ( MOD(Job,10)==0 ) THEN
      !
      !              COMPUTE INVERSE OF LOWER TRIANGULAR
      !
      DO kb = 1 , N
        k = N + 1 - kb
        Info = k
        IF ( T(k,k)==0.0E0 ) GOTO 99999
        T(k,k) = 1.0E0/T(k,k)
        temp = -T(k,k)
        IF ( k/=N ) CALL SSCAL(N-k,temp,T(k+1,k),1)
        km1 = k - 1
        IF ( km1>=1 ) THEN
          DO j = 1 , km1
            temp = T(k,j)
            T(k,j) = 0.0E0
            CALL SAXPY(N-k+1,temp,T(k,k),1,T(k,j),1)
          ENDDO
        ENDIF
      ENDDO
      Info = 0
    ELSE
      DO k = 1 , N
        Info = k
        IF ( T(k,k)==0.0E0 ) GOTO 99999
        T(k,k) = 1.0E0/T(k,k)
        temp = -T(k,k)
        CALL SSCAL(k-1,temp,T(1,k),1)
        kp1 = k + 1
        IF ( N>=kp1 ) THEN
          DO j = kp1 , N
            temp = T(k,j)
            T(k,j) = 0.0E0
            CALL SAXPY(k,temp,T(1,k),1,T(1,j),1)
          ENDDO
        ENDIF
      ENDDO
      Info = 0
    ENDIF
  ENDIF
  99999 END SUBROUTINE STRDI
