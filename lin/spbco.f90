!*==SPBCO.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SPBCO
      SUBROUTINE SPBCO(Abd,Lda,N,M,Rcond,Z,Info)
      IMPLICIT NONE
!*--SPBCO5
!***BEGIN PROLOGUE  SPBCO
!***PURPOSE  Factor a real symmetric positive definite matrix stored in
!            band form and estimate the condition number of the matrix.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B2
!***TYPE      SINGLE PRECISION (SPBCO-S, DPBCO-D, CPBCO-C)
!***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPBCO factors a real symmetric positive definite matrix
!     stored in band form and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, SPBFA is slightly faster.
!     To solve  A*X = B , follow SPBCO by SPBSL.
!     To compute  INVERSE(A)*C , follow SPBCO by SPBSL.
!     To compute  DETERMINANT(A) , follow SPBCO by SPBDI.
!
!     On Entry
!
!        ABD     REAL(LDA, N)
!                the matrix to be factored.  The columns of the upper
!                triangle are stored in the columns of ABD and the
!                diagonals of the upper triangle are stored in the
!                rows of ABD .  See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. M + 1 .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        M       INTEGER
!                the number of diagonals above the main diagonal.
!                0 .LE. M .LT. N .
!
!     On Return
!
!        ABD     an upper triangular matrix  R , stored in band
!                form, so that  A = TRANS(R)*R .
!                If  INFO .NE. 0 , the factorization is not complete.
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
!                underflows.  If INFO .NE. 0 , RCOND is unchanged.
!
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is singular to working precision, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                If  INFO .NE. 0 , Z  is unchanged.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!     Band Storage
!
!           If  A  is a symmetric positive definite band matrix,
!           the following program segment will set up the input.
!
!                   M = (band width above diagonal)
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses  M + 1  rows of  A , except for the  M by M
!           upper left triangle, which is ignored.
!
!     Example:  If the original matrix is
!
!           11 12 13  0  0  0
!           12 22 23 24  0  0
!           13 23 33 34 35  0
!            0 24 34 44 45 46
!            0  0 35 45 55 56
!            0  0  0 46 56 66
!
!     then  N = 6 , M = 2  and  ABD  should contain
!
!            *  * 13 24 35 46
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SASUM, SAXPY, SDOT, SPBFA, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPBCO
      INTEGER Lda , N , M , Info
      REAL Abd(Lda,*) , Z(*)
      REAL Rcond
!
      REAL SDOT , ek , t , wk , wkm
      REAL anorm , s , SASUM , sm , ynorm
      INTEGER i , j , j2 , k , kb , kp1 , l , la , lb , lm , mu
!
!     FIND NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  SPBCO
      DO j = 1 , N
        l = MIN(j,M+1)
        mu = MAX(M+2-j,1)
        Z(j) = SASUM(l,Abd(mu,j),1)
        k = j - l
        IF ( M>=mu ) THEN
          DO i = mu , M
            k = k + 1
            Z(k) = Z(k) + ABS(Abd(i,j))
          ENDDO
        ENDIF
      ENDDO
      anorm = 0.0E0
      DO j = 1 , N
        anorm = MAX(anorm,Z(j))
      ENDDO
!
!     FACTOR
!
      CALL SPBFA(Abd,Lda,N,M,Info)
      IF ( Info==0 ) THEN
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
        ek = 1.0E0
        DO j = 1 , N
          Z(j) = 0.0E0
        ENDDO
        DO k = 1 , N
          IF ( Z(k)/=0.0E0 ) ek = SIGN(ek,-Z(k))
          IF ( ABS(ek-Z(k))>Abd(M+1,k) ) THEN
            s = Abd(M+1,k)/ABS(ek-Z(k))
            CALL SSCAL(N,s,Z,1)
            ek = s*ek
          ENDIF
          wk = ek - Z(k)
          wkm = -ek - Z(k)
          s = ABS(wk)
          sm = ABS(wkm)
          wk = wk/Abd(M+1,k)
          wkm = wkm/Abd(M+1,k)
          kp1 = k + 1
          j2 = MIN(k+M,N)
          i = M + 1
          IF ( kp1<=j2 ) THEN
            DO j = kp1 , j2
              i = i - 1
              sm = sm + ABS(Z(j)+wkm*Abd(i,j))
              Z(j) = Z(j) + wk*Abd(i,j)
              s = s + ABS(Z(j))
            ENDDO
            IF ( s<sm ) THEN
              t = wkm - wk
              wk = wkm
              i = M + 1
              DO j = kp1 , j2
                i = i - 1
                Z(j) = Z(j) + t*Abd(i,j)
              ENDDO
            ENDIF
          ENDIF
          Z(k) = wk
        ENDDO
        s = 1.0E0/SASUM(N,Z,1)
        CALL SSCAL(N,s,Z,1)
!
!        SOLVE  R*Y = W
!
        DO kb = 1 , N
          k = N + 1 - kb
          IF ( ABS(Z(k))>Abd(M+1,k) ) THEN
            s = Abd(M+1,k)/ABS(Z(k))
            CALL SSCAL(N,s,Z,1)
          ENDIF
          Z(k) = Z(k)/Abd(M+1,k)
          lm = MIN(k-1,M)
          la = M + 1 - lm
          lb = k - lm
          t = -Z(k)
          CALL SAXPY(lm,t,Abd(la,k),1,Z(lb),1)
        ENDDO
        s = 1.0E0/SASUM(N,Z,1)
        CALL SSCAL(N,s,Z,1)
!
        ynorm = 1.0E0
!
!        SOLVE TRANS(R)*V = Y
!
        DO k = 1 , N
          lm = MIN(k-1,M)
          la = M + 1 - lm
          lb = k - lm
          Z(k) = Z(k) - SDOT(lm,Abd(la,k),1,Z(lb),1)
          IF ( ABS(Z(k))>Abd(M+1,k) ) THEN
            s = Abd(M+1,k)/ABS(Z(k))
            CALL SSCAL(N,s,Z,1)
            ynorm = s*ynorm
          ENDIF
          Z(k) = Z(k)/Abd(M+1,k)
        ENDDO
        s = 1.0E0/SASUM(N,Z,1)
        CALL SSCAL(N,s,Z,1)
        ynorm = s*ynorm
!
!        SOLVE  R*Z = W
!
        DO kb = 1 , N
          k = N + 1 - kb
          IF ( ABS(Z(k))>Abd(M+1,k) ) THEN
            s = Abd(M+1,k)/ABS(Z(k))
            CALL SSCAL(N,s,Z,1)
            ynorm = s*ynorm
          ENDIF
          Z(k) = Z(k)/Abd(M+1,k)
          lm = MIN(k-1,M)
          la = M + 1 - lm
          lb = k - lm
          t = -Z(k)
          CALL SAXPY(lm,t,Abd(la,k),1,Z(lb),1)
        ENDDO
!        MAKE ZNORM = 1.0
        s = 1.0E0/SASUM(N,Z,1)
        CALL SSCAL(N,s,Z,1)
        ynorm = s*ynorm
!
        IF ( anorm/=0.0E0 ) Rcond = ynorm/anorm
        IF ( anorm==0.0E0 ) Rcond = 0.0E0
      ENDIF
      END SUBROUTINE SPBCO
