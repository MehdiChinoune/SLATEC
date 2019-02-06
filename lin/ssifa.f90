!*==SSIFA.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SSIFA
      SUBROUTINE SSIFA(A,Lda,N,Kpvt,Info)
      IMPLICIT NONE
!*--SSIFA5
!***BEGIN PROLOGUE  SSIFA
!***PURPOSE  Factor a real symmetric matrix by elimination with
!            symmetric pivoting.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      SINGLE PRECISION (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     SSIFA factors a real symmetric matrix by elimination
!     with symmetric pivoting.
!
!     To solve  A*X = B , follow SSIFA by SSISL.
!     To compute  INVERSE(A)*C , follow SSIFA by SSISL.
!     To compute  DETERMINANT(A) , follow SSIFA by SSIDI.
!     To compute  INERTIA(A) , follow SSIFA by SSIDI.
!     To compute  INVERSE(A) , follow SSIFA by SSIDI.
!
!     On Entry
!
!        A       REAL(LDA,N)
!                the symmetric matrix to be factored.
!                Only the diagonal and upper triangle are used.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       a block diagonal matrix and the multipliers which
!                were used to obtain it.
!                The factorization can be written  A = U*D*TRANS(U)
!                where  U  is a product of permutation and unit
!                upper triangular matrices , TRANS(U) is the
!                transpose of  U , and  D  is block diagonal
!                with 1 by 1 and 2 by 2 blocks.
!
!        KPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if the K-th pivot block is singular.  This is
!                     not an error condition for this subroutine,
!                     but it does indicate that SSISL or SSIDI may
!                     divide by zero if called.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  ISAMAX, SAXPY, SSWAP
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
!***END PROLOGUE  SSIFA
      INTEGER Lda , N , Kpvt(*) , Info
      REAL A(Lda,*)
!
      REAL ak , akm1 , bk , bkm1 , denom , mulk , mulkm1 , t
      REAL absakk , alpha , colmax , rowmax
      INTEGER imax , imaxp1 , j , jj , jmax , k , km1 , km2 , kstep , ISAMAX
      LOGICAL swap
!***FIRST EXECUTABLE STATEMENT  SSIFA
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
!
      alpha = (1.0E0+SQRT(17.0E0))/8.0E0
!
      Info = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
      k = N
!
!        LEAVE THE LOOP IF K=0 OR K=1.
!
      DO WHILE ( k/=0 )
        IF ( k>1 ) THEN
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
!        REQUIRED.
!
          km1 = k - 1
          absakk = ABS(A(k,k))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
          imax = ISAMAX(k-1,A(1,k),1)
          colmax = ABS(A(imax,k))
          IF ( absakk<alpha*colmax ) THEN
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
            rowmax = 0.0E0
            imaxp1 = imax + 1
            DO j = imaxp1 , k
              rowmax = MAX(rowmax,ABS(A(imax,j)))
            ENDDO
            IF ( imax/=1 ) THEN
              jmax = ISAMAX(imax-1,A(1,imax),1)
              rowmax = MAX(rowmax,ABS(A(jmax,imax)))
            ENDIF
            IF ( ABS(A(imax,imax))>=alpha*rowmax ) THEN
              kstep = 1
              swap = .TRUE.
            ELSEIF ( absakk<alpha*colmax*(colmax/rowmax) ) THEN
              kstep = 2
              swap = imax/=km1
            ELSE
              kstep = 1
              swap = .FALSE.
            ENDIF
          ELSE
            kstep = 1
            swap = .FALSE.
          ENDIF
          IF ( MAX(absakk,colmax)==0.0E0 ) THEN
!
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
            Kpvt(k) = k
            Info = k
          ELSEIF ( kstep==2 ) THEN
!
!           2 X 2 PIVOT BLOCK.
!
            IF ( swap ) THEN
!
!              PERFORM AN INTERCHANGE.
!
              CALL SSWAP(imax,A(1,imax),1,A(1,k-1),1)
              DO jj = imax , km1
                j = km1 + imax - jj
                t = A(j,k-1)
                A(j,k-1) = A(imax,j)
                A(imax,j) = t
              ENDDO
              t = A(k-1,k)
              A(k-1,k) = A(imax,k)
              A(imax,k) = t
            ENDIF
!
!           PERFORM THE ELIMINATION.
!
            km2 = k - 2
            IF ( km2/=0 ) THEN
              ak = A(k,k)/A(k-1,k)
              akm1 = A(k-1,k-1)/A(k-1,k)
              denom = 1.0E0 - ak*akm1
              DO jj = 1 , km2
                j = km1 - jj
                bk = A(j,k)/A(k-1,k)
                bkm1 = A(j,k-1)/A(k-1,k)
                mulk = (akm1*bk-bkm1)/denom
                mulkm1 = (ak*bkm1-bk)/denom
                t = mulk
                CALL SAXPY(j,t,A(1,k),1,A(1,j),1)
                t = mulkm1
                CALL SAXPY(j,t,A(1,k-1),1,A(1,j),1)
                A(j,k) = mulk
                A(j,k-1) = mulkm1
              ENDDO
            ENDIF
!
!           SET THE PIVOT ARRAY.
!
            Kpvt(k) = 1 - k
            IF ( swap ) Kpvt(k) = -imax
            Kpvt(k-1) = Kpvt(k)
          ELSE
!
!           1 X 1 PIVOT BLOCK.
!
            IF ( swap ) THEN
!
!              PERFORM AN INTERCHANGE.
!
              CALL SSWAP(imax,A(1,imax),1,A(1,k),1)
              DO jj = imax , k
                j = k + imax - jj
                t = A(j,k)
                A(j,k) = A(imax,j)
                A(imax,j) = t
              ENDDO
            ENDIF
!
!           PERFORM THE ELIMINATION.
!
            DO jj = 1 , km1
              j = k - jj
              mulk = -A(j,k)/A(k,k)
              t = mulk
              CALL SAXPY(j,t,A(1,k),1,A(1,j),1)
              A(j,k) = mulk
            ENDDO
!
!           SET THE PIVOT ARRAY.
!
            Kpvt(k) = k
            IF ( swap ) Kpvt(k) = imax
          ENDIF
          k = k - kstep
        ELSE
          Kpvt(1) = 1
          IF ( A(1,1)==0.0E0 ) Info = 1
          EXIT
        ENDIF
      ENDDO
      END SUBROUTINE SSIFA
