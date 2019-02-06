!*==CGECO.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CGECO
      SUBROUTINE CGECO(A,Lda,N,Ipvt,Rcond,Z)
      IMPLICIT NONE
!*--CGECO5
!***BEGIN PROLOGUE  CGECO
!***PURPOSE  Factor a matrix using Gaussian elimination and estimate
!            the condition number of the matrix.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C1
!***TYPE      COMPLEX (SGECO-S, DGECO-D, CGECO-C)
!***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CGECO factors a complex matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, CGEFA is slightly faster.
!     To solve  A*X = B , follow CGECO By CGESL.
!     To Compute  INVERSE(A)*C , follow CGECO by CGESL.
!     To compute  DETERMINANT(A) , follow CGECO by CGEDI.
!     To compute  INVERSE(A) , follow CGECO by CGEDI.
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       COMPLEX(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC, CGEFA, CSSCAL, SCASUM
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CGECO
      INTEGER Lda , N , Ipvt(*)
      COMPLEX A(Lda,*) , Z(*)
      REAL Rcond
!
      COMPLEX CDOTC , ek , t , wk , wkm
      REAL anorm , s , SCASUM , sm , ynorm
      INTEGER info , j , k , kb , kp1 , l
      COMPLEX zdum , zdum1 , zdum2 , CSIGN1
      REAL CABS1
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
      CSIGN1(zdum1,zdum2) = CABS1(zdum1)*(zdum2/CABS1(zdum2))
!
!     COMPUTE 1-NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  CGECO
      anorm = 0.0E0
      DO j = 1 , N
        anorm = MAX(anorm,SCASUM(N,A(1,j),1))
      ENDDO
!
!     FACTOR
!
      CALL CGEFA(A,Lda,N,Ipvt,info)
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
!     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE CTRANS(U)*W = E
!
      ek = (1.0E0,0.0E0)
      DO j = 1 , N
        Z(j) = (0.0E0,0.0E0)
      ENDDO
      DO k = 1 , N
        IF ( CABS1(Z(k))/=0.0E0 ) ek = CSIGN1(ek,-Z(k))
        IF ( CABS1(ek-Z(k))>CABS1(A(k,k)) ) THEN
          s = CABS1(A(k,k))/CABS1(ek-Z(k))
          CALL CSSCAL(N,s,Z,1)
          ek = CMPLX(s,0.0E0)*ek
        ENDIF
        wk = ek - Z(k)
        wkm = -ek - Z(k)
        s = CABS1(wk)
        sm = CABS1(wkm)
        IF ( CABS1(A(k,k))==0.0E0 ) THEN
          wk = (1.0E0,0.0E0)
          wkm = (1.0E0,0.0E0)
        ELSE
          wk = wk/CONJG(A(k,k))
          wkm = wkm/CONJG(A(k,k))
        ENDIF
        kp1 = k + 1
        IF ( kp1<=N ) THEN
          DO j = kp1 , N
            sm = sm + CABS1(Z(j)+wkm*CONJG(A(k,j)))
            Z(j) = Z(j) + wk*CONJG(A(k,j))
            s = s + CABS1(Z(j))
          ENDDO
          IF ( s<sm ) THEN
            t = wkm - wk
            wk = wkm
            DO j = kp1 , N
              Z(j) = Z(j) + t*CONJG(A(k,j))
            ENDDO
          ENDIF
        ENDIF
        Z(k) = wk
      ENDDO
      s = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,s,Z,1)
!
!     SOLVE CTRANS(L)*Y = W
!
      DO kb = 1 , N
        k = N + 1 - kb
        IF ( k<N ) Z(k) = Z(k) + CDOTC(N-k,A(k+1,k),1,Z(k+1),1)
        IF ( CABS1(Z(k))>1.0E0 ) THEN
          s = 1.0E0/CABS1(Z(k))
          CALL CSSCAL(N,s,Z,1)
        ENDIF
        l = Ipvt(k)
        t = Z(l)
        Z(l) = Z(k)
        Z(k) = t
      ENDDO
      s = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,s,Z,1)
!
      ynorm = 1.0E0
!
!     SOLVE L*V = Y
!
      DO k = 1 , N
        l = Ipvt(k)
        t = Z(l)
        Z(l) = Z(k)
        Z(k) = t
        IF ( k<N ) CALL CAXPY(N-k,t,A(k+1,k),1,Z(k+1),1)
        IF ( CABS1(Z(k))>1.0E0 ) THEN
          s = 1.0E0/CABS1(Z(k))
          CALL CSSCAL(N,s,Z,1)
          ynorm = s*ynorm
        ENDIF
      ENDDO
      s = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,s,Z,1)
      ynorm = s*ynorm
!
!     SOLVE  U*Z = V
!
      DO kb = 1 , N
        k = N + 1 - kb
        IF ( CABS1(Z(k))>CABS1(A(k,k)) ) THEN
          s = CABS1(A(k,k))/CABS1(Z(k))
          CALL CSSCAL(N,s,Z,1)
          ynorm = s*ynorm
        ENDIF
        IF ( CABS1(A(k,k))/=0.0E0 ) Z(k) = Z(k)/A(k,k)
        IF ( CABS1(A(k,k))==0.0E0 ) Z(k) = (1.0E0,0.0E0)
        t = -Z(k)
        CALL CAXPY(k-1,t,A(1,k),1,Z(1),1)
      ENDDO
!     MAKE ZNORM = 1.0
      s = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,s,Z,1)
      ynorm = s*ynorm
!
      IF ( anorm/=0.0E0 ) Rcond = ynorm/anorm
      IF ( anorm==0.0E0 ) Rcond = 0.0E0
      END SUBROUTINE CGECO
