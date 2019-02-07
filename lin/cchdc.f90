!*==CCHDC.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CCHDC
SUBROUTINE CCHDC(A,Lda,P,Work,Jpvt,Job,Info)
  IMPLICIT NONE
  !*--CCHDC5
  !***BEGIN PROLOGUE  CCHDC
  !***PURPOSE  Compute the Cholesky decomposition of a positive definite
  !            matrix.  A pivoting option allows the user to estimate the
  !            condition number of a positive definite matrix or determine
  !            the rank of a positive semidefinite matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D1B
  !***TYPE      COMPLEX (SCHDC-S, DCHDC-D, CCHDC-C)
  !***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX,
  !             POSITIVE DEFINITE
  !***AUTHOR  Dongarra, J., (ANL)
  !           Stewart, G. W., (U. of Maryland)
  !***DESCRIPTION
  !
  !     CCHDC computes the Cholesky decomposition of a positive definite
  !     matrix.  A pivoting option allows the user to estimate the
  !     condition of a positive definite matrix or determine the rank
  !     of a positive semidefinite matrix.
  !
  !     On Entry
  !
  !         A      COMPLEX(LDA,P).
  !                A contains the matrix whose decomposition is to
  !                be computed.  Only the upper half of A need be stored.
  !                The lower part of The array A is not referenced.
  !
  !         LDA    INTEGER.
  !                LDA is the leading dimension of the array A.
  !
  !         P      INTEGER.
  !                P is the order of the matrix.
  !
  !         WORK   COMPLEX.
  !                WORK is a work array.
  !
  !         JPVT   INTEGER(P).
  !                JPVT contains integers that control the selection
  !                of the pivot elements, if pivoting has been requested.
  !                Each diagonal element A(K,K)
  !                is placed in one of three classes according to the
  !                value of JPVT(K)).
  !
  !                   If JPVT(K)) .GT. 0, then X(K) is an initial
  !                                      element.
  !
  !                   If JPVT(K)) .EQ. 0, then X(K) is a free element.
  !
  !                   If JPVT(K)) .LT. 0, then X(K) is a final element.
  !
  !                Before the decomposition is computed, initial elements
  !                are moved by symmetric row and column interchanges to
  !                the beginning of the array A and final
  !                elements to the end.  Both initial and final elements
  !                are frozen in place during the computation and only
  !                free elements are moved.  At the K-th stage of the
  !                reduction, if A(K,K) is occupied by a free element
  !                it is interchanged with the largest free element
  !                A(L,L) with L .GE. K.  JPVT is not referenced if
  !                JOB .EQ. 0.
  !
  !        JOB     INTEGER.
  !                JOB is an integer that initiates column pivoting.
  !                IF JOB .EQ. 0, no pivoting is done.
  !                IF JOB .NE. 0, pivoting is done.
  !
  !     On Return
  !
  !         A      A contains in its upper half the Cholesky factor
  !                of the matrix A as it has been permuted by pivoting.
  !
  !         JPVT   JPVT(J) contains the index of the diagonal element
  !                of A that was moved into the J-th position,
  !                provided pivoting was requested.
  !
  !         INFO   contains the index of the last positive diagonal
  !                element of the Cholesky factor.
  !
  !     For positive definite matrices INFO = P is the normal return.
  !     For pivoting with positive semidefinite matrices INFO will
  !     in general be less than P.  However, INFO may be greater than
  !     the rank of A, since rounding error can cause an otherwise zero
  !     element to be positive.  Indefinite systems will always cause
  !     INFO to be less than P.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CAXPY, CSWAP
  !***REVISION HISTORY  (YYMMDD)
  !   790319  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CCHDC
  INTEGER Lda, P, Jpvt(*), Job, Info
  COMPLEX A(Lda,*), Work(*)
  !
  INTEGER pu, pl, plp1, j, jp, jt, k, kb, km1, kp1, l, maxl
  COMPLEX temp
  REAL maxdia
  LOGICAL swapk, negk
  !***FIRST EXECUTABLE STATEMENT  CCHDC
  pl = 1
  pu = 0
  Info = P
  IF ( Job/=0 ) THEN
    !
    !        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
    !        THE ELEMENTS ACCORDING TO JPVT.
    !
    DO k = 1, P
      swapk = Jpvt(k)>0
      negk = Jpvt(k)<0
      Jpvt(k) = k
      IF ( negk ) Jpvt(k) = -Jpvt(k)
      IF ( swapk ) THEN
        IF ( k/=pl ) THEN
          CALL CSWAP(pl-1,A(1,k),1,A(1,pl),1)
          temp = A(k,k)
          A(k,k) = A(pl,pl)
          A(pl,pl) = temp
          A(pl,k) = CONJG(A(pl,k))
          plp1 = pl + 1
          IF ( P>=plp1 ) THEN
            DO j = plp1, P
              IF ( j<k ) THEN
                temp = CONJG(A(pl,j))
                A(pl,j) = CONJG(A(j,k))
                A(j,k) = temp
              ELSEIF ( j/=k ) THEN
                temp = A(k,j)
                A(k,j) = A(pl,j)
                A(pl,j) = temp
              ENDIF
            ENDDO
          ENDIF
          Jpvt(k) = Jpvt(pl)
          Jpvt(pl) = k
        ENDIF
        pl = pl + 1
      ENDIF
    ENDDO
    pu = P
    IF ( P>=pl ) THEN
      DO kb = pl, P
        k = P - kb + pl
        IF ( Jpvt(k)<0 ) THEN
          Jpvt(k) = -Jpvt(k)
          IF ( pu/=k ) THEN
            CALL CSWAP(k-1,A(1,k),1,A(1,pu),1)
            temp = A(k,k)
            A(k,k) = A(pu,pu)
            A(pu,pu) = temp
            A(k,pu) = CONJG(A(k,pu))
            kp1 = k + 1
            IF ( P>=kp1 ) THEN
              DO j = kp1, P
                IF ( j<pu ) THEN
                  temp = CONJG(A(k,j))
                  A(k,j) = CONJG(A(j,pu))
                  A(j,pu) = temp
                ELSEIF ( j/=pu ) THEN
                  temp = A(k,j)
                  A(k,j) = A(pu,j)
                  A(pu,j) = temp
                ENDIF
              ENDDO
            ENDIF
            jt = Jpvt(k)
            Jpvt(k) = Jpvt(pu)
            Jpvt(pu) = jt
          ENDIF
          pu = pu - 1
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  DO k = 1, P
    !
    !        REDUCTION LOOP.
    !
    maxdia = REAL(A(k,k))
    kp1 = k + 1
    maxl = k
    !
    !        DETERMINE THE PIVOT ELEMENT.
    !
    IF ( k>=pl.AND.k<pu ) THEN
      DO l = kp1, pu
        IF ( REAL(A(l,l))>maxdia ) THEN
          maxdia = REAL(A(l,l))
          maxl = l
        ENDIF
      ENDDO
    ENDIF
    !
    !        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.
    !
    IF ( maxdia>0.0E0 ) THEN
      IF ( k/=maxl ) THEN
        !
        !           START THE PIVOTING AND UPDATE JPVT.
        !
        km1 = k - 1
        CALL CSWAP(km1,A(1,k),1,A(1,maxl),1)
        A(maxl,maxl) = A(k,k)
        A(k,k) = CMPLX(maxdia,0.0E0)
        jp = Jpvt(maxl)
        Jpvt(maxl) = Jpvt(k)
        Jpvt(k) = jp
        A(k,maxl) = CONJG(A(k,maxl))
      ENDIF
      !
      !        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
      !
      Work(k) = CMPLX(SQRT(REAL(A(k,k))),0.0E0)
      A(k,k) = Work(k)
      IF ( P>=kp1 ) THEN
        DO j = kp1, P
          IF ( k/=maxl ) THEN
            IF ( j<maxl ) THEN
              temp = CONJG(A(k,j))
              A(k,j) = CONJG(A(j,maxl))
              A(j,maxl) = temp
            ELSEIF ( j/=maxl ) THEN
              temp = A(k,j)
              A(k,j) = A(maxl,j)
              A(maxl,j) = temp
            ENDIF
          ENDIF
          A(k,j) = A(k,j)/Work(k)
          Work(j) = CONJG(A(k,j))
          temp = -A(k,j)
          CALL CAXPY(j-k,temp,Work(kp1),1,A(kp1,j),1)
        ENDDO
      ENDIF
    ELSE
      Info = k - 1
      EXIT
    ENDIF
  ENDDO
END SUBROUTINE CCHDC
