!** DSIFA
SUBROUTINE DSIFA(A,Lda,N,Kpvt,Info)
  !>
  !***
  !  Factor a real symmetric matrix by elimination with
  !            symmetric pivoting.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1A
  !***
  ! **Type:**      DOUBLE PRECISION (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, SYMMETRIC
  !***
  ! **Author:**  Bunch, J., (UCSD)
  !***
  ! **Description:**
  !
  !     DSIFA factors a double precision symmetric matrix by elimination
  !     with symmetric pivoting.
  !
  !     To solve  A*X = B, follow DSIFA by DSISL.
  !     To compute  INVERSE(A)*C, follow DSIFA by DSISL.
  !     To compute  DETERMINANT(A), follow DSIFA by DSIDI.
  !     To compute  INERTIA(A), follow DSIFA by DSIDI.
  !     To compute  INVERSE(A), follow DSIFA by DSIDI.
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA,N)
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
  !                upper triangular matrices, TRANS(U) is the
  !                transpose of  U, and  D  is block diagonal
  !                with 1 by 1 and 2 by 2 blocks.
  !
  !        KPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if the K-th pivot block is singular.  This is
  !                     not an error condition for this subroutine,
  !                     but it does indicate that DSISL or DSIDI may
  !                     divide by zero if called.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DSWAP, IDAMAX

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Modified routine equivalence list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Lda, N, Kpvt(*), Info
  REAL(8) :: A(Lda,*)
  !
  REAL(8) :: ak, akm1, bk, bkm1, denom, mulk, mulkm1, t
  REAL(8) :: absakk, alpha, colmax, rowmax
  INTEGER imax, imaxp1, j, jj, jmax, k, km1, km2, kstep
  LOGICAL swap
  !* FIRST EXECUTABLE STATEMENT  DSIFA
  !
  !     INITIALIZE
  !
  !     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
  !
  alpha = (1.0D0+SQRT(17.0D0))/8.0D0
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
      imax = IDAMAX(k-1,A(1,k),1)
      colmax = ABS(A(imax,k))
      IF ( absakk<alpha*colmax ) THEN
        !
        !           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
        !           ROW IMAX.
        !
        rowmax = 0.0D0
        imaxp1 = imax + 1
        DO j = imaxp1, k
          rowmax = MAX(rowmax,ABS(A(imax,j)))
        END DO
        IF ( imax/=1 ) THEN
          jmax = IDAMAX(imax-1,A(1,imax),1)
          rowmax = MAX(rowmax,ABS(A(jmax,imax)))
        END IF
        IF ( ABS(A(imax,imax))>=alpha*rowmax ) THEN
          kstep = 1
          swap = .TRUE.
        ELSEIF ( absakk<alpha*colmax*(colmax/rowmax) ) THEN
          kstep = 2
          swap = imax/=km1
        ELSE
          kstep = 1
          swap = .FALSE.
        END IF
      ELSE
        kstep = 1
        swap = .FALSE.
      END IF
      IF ( MAX(absakk,colmax)==0.0D0 ) THEN
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
          CALL DSWAP(imax,A(1,imax),1,A(1,k-1),1)
          DO jj = imax, km1
            j = km1 + imax - jj
            t = A(j,k-1)
            A(j,k-1) = A(imax,j)
            A(imax,j) = t
          END DO
          t = A(k-1,k)
          A(k-1,k) = A(imax,k)
          A(imax,k) = t
        END IF
        !
        !           PERFORM THE ELIMINATION.
        !
        km2 = k - 2
        IF ( km2/=0 ) THEN
          ak = A(k,k)/A(k-1,k)
          akm1 = A(k-1,k-1)/A(k-1,k)
          denom = 1.0D0 - ak*akm1
          DO jj = 1, km2
            j = km1 - jj
            bk = A(j,k)/A(k-1,k)
            bkm1 = A(j,k-1)/A(k-1,k)
            mulk = (akm1*bk-bkm1)/denom
            mulkm1 = (ak*bkm1-bk)/denom
            t = mulk
            CALL DAXPY(j,t,A(1,k),1,A(1,j),1)
            t = mulkm1
            CALL DAXPY(j,t,A(1,k-1),1,A(1,j),1)
            A(j,k) = mulk
            A(j,k-1) = mulkm1
          END DO
        END IF
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
          CALL DSWAP(imax,A(1,imax),1,A(1,k),1)
          DO jj = imax, k
            j = k + imax - jj
            t = A(j,k)
            A(j,k) = A(imax,j)
            A(imax,j) = t
          END DO
        END IF
        !
        !           PERFORM THE ELIMINATION.
        !
        DO jj = 1, km1
          j = k - jj
          mulk = -A(j,k)/A(k,k)
          t = mulk
          CALL DAXPY(j,t,A(1,k),1,A(1,j),1)
          A(j,k) = mulk
        END DO
        !
        !           SET THE PIVOT ARRAY.
        !
        Kpvt(k) = k
        IF ( swap ) Kpvt(k) = imax
      END IF
      k = k - kstep
    ELSE
      Kpvt(1) = 1
      IF ( A(1,1)==0.0D0 ) Info = 1
      EXIT
    END IF
  END DO
END SUBROUTINE DSIFA
