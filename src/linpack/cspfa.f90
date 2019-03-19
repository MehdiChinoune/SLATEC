!** CSPFA
SUBROUTINE CSPFA(Ap,N,Kpvt,Info)
  IMPLICIT NONE
  !>
  !***
  !  Factor a complex symmetric matrix stored in packed form by
  !            elimination with symmetric pivoting.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2C1
  !***
  ! **Type:**      COMPLEX (SSPFA-S, DSPFA-D, CHPFA-C, CSPFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
  !             SYMMETRIC
  !***
  ! **Author:**  Bunch, J., (UCSD)
  !***
  ! **Description:**
  !
  !     CSPFA factors a complex symmetric matrix stored in
  !     packed form by elimination with symmetric pivoting.
  !
  !     To solve  A*X = B, follow CSPFA by CSPSL.
  !     To compute  INVERSE(A)*C, follow CSPFA by CSPSL.
  !     To compute  DETERMINANT(A), follow CSPFA by CSPDI.
  !     To compute  INVERSE(A), follow CSPFA by CSPDI.
  !
  !     On Entry
  !
  !        AP      COMPLEX (N*(N+1)/2)
  !                the packed form of a symmetric matrix  A .  The
  !                columns of the upper triangle are stored sequentially
  !                in a one-dimensional array of length  N*(N+1)/2 .
  !                See comments below for details.
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        AP      a block diagonal matrix and the multipliers which
  !                were used to obtain it stored in packed form.
  !                The factorization can be written  A = U*D*TRANS(U)
  !                where  U  is a product of permutation and unit
  !                upper triangular matrices, TRANS(U) is the
  !                transpose of  U, and  D  is block diagonal
  !                with 1 by 1 and 2 by 2 blocks.
  !
  !        KVPT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if the K-th pivot block is singular.  This is
  !                     not an error condition for this subroutine,
  !                     but it does indicate that CSPSL or CSPDI may
  !                     divide by zero if called.
  !
  !     Packed Storage
  !
  !          The following program segment will pack the upper
  !          triangle of a symmetric matrix.
  !
  !                K = 0
  !                DO 20 J = 1, N
  !                   DO 10 I = 1, J
  !                      K = K + 1
  !                      AP(K)  = A(I,J)
  !             10    CONTINUE
  !             20 CONTINUE
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CAXPY, CSWAP, ICAMAX

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891107  Corrected category and modified routine equivalence
  !           list.  (WRB)
  !   891107  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER N, Kpvt(*), Info
  COMPLEX Ap(*)
  !
  COMPLEX ak, akm1, bk, bkm1, denom, mulk, mulkm1, t
  REAL absakk, alpha, colmax, rowmax
  INTEGER ICAMAX, ij, ik, ikm1, im, imax, imaxp1, imim, imj, imk
  INTEGER j, jj, jk, jkm1, jmax, jmim, k, kk, km1, km1k, km1km1, &
    km2, kstep
  LOGICAL swap
  REAL, EXTERNAL :: CABS1
  !* FIRST EXECUTABLE STATEMENT  CSPFA
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
  ik = (N*(N-1))/2
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
      kk = ik + k
      absakk = CABS1(Ap(kk))
      !
      !        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
      !        COLUMN K.
      !
      imax = ICAMAX(k-1,Ap(ik+1),1)
      imk = ik + imax
      colmax = CABS1(Ap(imk))
      IF ( absakk<alpha*colmax ) THEN
        !
        !           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
        !           ROW IMAX.
        !
        rowmax = 0.0E0
        imaxp1 = imax + 1
        im = imax*(imax-1)/2
        imj = im + 2*imax
        DO j = imaxp1, k
          rowmax = MAX(rowmax,CABS1(Ap(imj)))
          imj = imj + j
        ENDDO
        IF ( imax/=1 ) THEN
          jmax = ICAMAX(imax-1,Ap(im+1),1)
          jmim = jmax + im
          rowmax = MAX(rowmax,CABS1(Ap(jmim)))
        ENDIF
        imim = imax + im
        IF ( CABS1(Ap(imim))>=alpha*rowmax ) THEN
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
        km1k = ik + k - 1
        ikm1 = ik - (k-1)
        IF ( swap ) THEN
          !
          !              PERFORM AN INTERCHANGE.
          !
          CALL CSWAP(imax,Ap(im+1),1,Ap(ikm1+1),1)
          imj = ikm1 + imax
          DO jj = imax, km1
            j = km1 + imax - jj
            jkm1 = ikm1 + j
            t = Ap(jkm1)
            Ap(jkm1) = Ap(imj)
            Ap(imj) = t
            imj = imj - (j-1)
          ENDDO
          t = Ap(km1k)
          Ap(km1k) = Ap(imk)
          Ap(imk) = t
        ENDIF
        !
        !           PERFORM THE ELIMINATION.
        !
        km2 = k - 2
        IF ( km2/=0 ) THEN
          ak = Ap(kk)/Ap(km1k)
          km1km1 = ikm1 + k - 1
          akm1 = Ap(km1km1)/Ap(km1k)
          denom = 1.0E0 - ak*akm1
          ij = ik - (k-1) - (k-2)
          DO jj = 1, km2
            j = km1 - jj
            jk = ik + j
            bk = Ap(jk)/Ap(km1k)
            jkm1 = ikm1 + j
            bkm1 = Ap(jkm1)/Ap(km1k)
            mulk = (akm1*bk-bkm1)/denom
            mulkm1 = (ak*bkm1-bk)/denom
            t = mulk
            CALL CAXPY(j,t,Ap(ik+1),1,Ap(ij+1),1)
            t = mulkm1
            CALL CAXPY(j,t,Ap(ikm1+1),1,Ap(ij+1),1)
            Ap(jk) = mulk
            Ap(jkm1) = mulkm1
            ij = ij - (j-1)
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
          CALL CSWAP(imax,Ap(im+1),1,Ap(ik+1),1)
          imj = ik + imax
          DO jj = imax, k
            j = k + imax - jj
            jk = ik + j
            t = Ap(jk)
            Ap(jk) = Ap(imj)
            Ap(imj) = t
            imj = imj - (j-1)
          ENDDO
        ENDIF
        !
        !           PERFORM THE ELIMINATION.
        !
        ij = ik - (k-1)
        DO jj = 1, km1
          j = k - jj
          jk = ik + j
          mulk = -Ap(jk)/Ap(kk)
          t = mulk
          CALL CAXPY(j,t,Ap(ik+1),1,Ap(ij+1),1)
          Ap(jk) = mulk
          ij = ij - (j-1)
        ENDDO
        !
        !           SET THE PIVOT ARRAY.
        !
        Kpvt(k) = k
        IF ( swap ) Kpvt(k) = imax
      ENDIF
      ik = ik - (k-1)
      IF ( kstep==2 ) ik = ik - (k-2)
      k = k - kstep
    ELSE
      Kpvt(1) = 1
      IF ( CABS1(Ap(1))==0.0E0 ) Info = 1
      EXIT
    ENDIF
  ENDDO
END SUBROUTINE CSPFA
