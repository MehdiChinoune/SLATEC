!** DMGSBV
PURE SUBROUTINE DMGSBV(M,N,A,Ia,Niv,Iflag,S,P,Ip,Inhomo,V,W,Wcnd)
  !> Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (MGSBV-S, DMGSBV-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  ! Orthogonalize a set of N double precision vectors and determine their
  ! rank.
  !
  !- *********************************************************************
  ! INPUT
  !- *********************************************************************
  !   M = dimension of vectors.
  !   N = no. of vectors.
  !   A = array whose first N cols contain the vectors.
  !   IA = first dimension of array A (col length).
  !   NIV = number of independent vectors needed.
  !   INHOMO = 1 corresponds to having a non-zero particular solution.
  !   V = particular solution vector (not included in the pivoting).
  !   INDPVT = 1 means pivoting will not be used.
  !
  !- *********************************************************************
  ! OUTPUT
  !- *********************************************************************
  !   NIV = no. of linear independent vectors in input set.
  !     A = matrix whose first NIV cols. contain NIV orthogonal vectors
  !         which span the vector space determined by the input vectors.
  !   IFLAG
  !          = 0 success
  !          = 1 incorrect input
  !          = 2 rank of new vectors less than N
  !   P = decomposition matrix.  P is upper triangular and
  !             (old vectors) = (new vectors) * P.
  !         The old vectors will be reordered due to pivoting.
  !         The dimension of P must be >= N*(N+1)/2.
  !             (  N*(2*N+1) when N /= NFCC )
  !   IP = pivoting vector. The dimension of IP must be >= N.
  !             (  2*N when N /= NFCC )
  !   S = square of norms of incoming vectors.
  !   V = vector which is orthogonal to the vectors of A.
  !   W = orthogonalization information for the vector V.
  !   WCND = worst case (smallest) norm decrement value of the
  !          vectors being orthogonalized  (represents a test
  !          for linear dependence of the vectors).
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DDOT, DPRVEC
  !***
  ! COMMON BLOCKS    DML18J, DML5MC

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   890921  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : sru_com, eps_com, indpvt_com, nfcc_com
  !
  INTEGER, INTENT(IN) :: Ia, Inhomo, M, N
  INTEGER, INTENT(INOUT) :: Niv
  INTEGER, INTENT(OUT) :: Iflag, Ip(2*N)
  REAL(DP), INTENT(INOUT) :: A(Ia,M), V(M)
  REAL(DP), INTENT(OUT) :: Wcnd, P(N*(2*N+1)), S(M), W(M)
  !
  INTEGER :: i, ip1, ix, iz, j, jk, jp, jq, jy, jz, k, kd, kj, kp, l, lix, lr, &
    m2, nivn, nmnr, nn, np1, nr, nrm1
  REAL(DP) :: dot, pjp, psave, ry, sv, t, vl, vnorm, y
  !
  !* FIRST EXECUTABLE STATEMENT  DMGSBV
  IF( M>0 .AND. N>0 .AND. Ia>=M ) THEN
    !        BEGIN BLOCK PERMITTING ...EXITS TO 270
    !           BEGIN BLOCK PERMITTING ...EXITS TO 260
    !
    jp = 0
    Iflag = 0
    np1 = N + 1
    y = 0._DP
    m2 = M/2
    !
    !              CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH
    !              FOR VECTOR WITH LARGEST MAGNITUDE
    !
    j = 0
    DO i = 1, N
      vl = NORM2(A(1:M,i))**2
      S(i) = vl
      IF( N/=nfcc_com ) THEN
        j = 2*i - 1
        P(j) = vl
        Ip(j) = j
      END IF
      j = j + 1
      P(j) = vl
      Ip(j) = j
      IF( vl>y ) THEN
        y = vl
        ix = i
      END IF
    END DO
    IF( indpvt_com==1 ) THEN
      ix = 1
      y = P(1)
    END IF
    lix = ix
    IF( N/=nfcc_com ) lix = 2*ix - 1
    P(lix) = P(1)
    S(np1) = 0._DP
    IF( Inhomo==1 ) S(np1) = NORM2(V(1:M))**2
    Wcnd = 1._DP
    nivn = Niv
    Niv = 0
    !
    !           ...EXIT
    IF( y/=0._DP ) THEN
      !              *********************************************************
      DO nr = 1, N
        !                 BEGIN BLOCK PERMITTING ...EXITS TO 230
        !              ......EXIT
        IF( nivn==Niv ) EXIT
        Niv = nr
        IF( ix/=nr ) THEN
          !
          !                       PIVOTING OF COLUMNS OF P MATRIX
          !
          nn = N
          lix = ix
          lr = nr
          IF( N/=nfcc_com ) THEN
            nn = nfcc_com
            lix = 2*ix - 1
            lr = 2*nr - 1
          END IF
          IF( nr/=1 ) THEN
            kd = lix - lr
            kj = lr
            nrm1 = lr - 1
            DO j = 1, nrm1
              psave = P(kj)
              jk = kj + kd
              P(kj) = P(jk)
              P(jk) = psave
              kj = kj + nn - j
            END DO
            jy = jk + nmnr
            jz = jy - kd
            P(jy) = P(jz)
          END IF
          iz = Ip(lix)
          Ip(lix) = Ip(lr)
          Ip(lr) = iz
          sv = S(ix)
          S(ix) = S(nr)
          S(nr) = sv
          IF( N/=nfcc_com ) THEN
            IF( nr/=1 ) THEN
              kj = lr + 1
              DO k = 1, nrm1
                psave = P(kj)
                jk = kj + kd
                P(kj) = P(jk)
                P(jk) = psave
                kj = kj + nfcc_com - k
              END DO
            END IF
            iz = Ip(lix+1)
            Ip(lix+1) = Ip(lr+1)
            Ip(lr+1) = iz
          END IF
          !
          !                       PIVOTING OF COLUMNS OF VECTORS
          !
          DO l = 1, M
            t = A(l,ix)
            A(l,ix) = A(l,nr)
            A(l,nr) = t
          END DO
        END IF
        !
        !                    CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL
        !                    VECTOR
        !
        jp = jp + 1
        P(jp) = y
        ry = 1._DP/y
        nmnr = N - nr
        IF( N/=nfcc_com ) THEN
          nmnr = nfcc_com - (2*nr-1)
          jp = jp + 1
          P(jp) = 0._DP
          kp = jp + nmnr
          P(kp) = y
        END IF
        IF( nr/=N .AND. nivn/=Niv ) THEN
          !
          !                       CALCULATE ORTHOGONAL PROJECTION VECTORS AND
          !                       SEARCH FOR LARGEST NORM
          !
          y = 0._DP
          ip1 = nr + 1
          ix = ip1
          !                       ************************************************
          DO j = ip1, N
            dot = DOT_PRODUCT(A(1:M,nr),A(1:M,j))
            jp = jp + 1
            jq = jp + nmnr
            IF( N/=nfcc_com ) jq = jq + nmnr - 1
            P(jq) = P(jp) - dot*(dot*ry)
            P(jp) = dot*ry
            DO i = 1, M
              A(i,j) = A(i,j) - P(jp)*A(i,nr)
            END DO
            IF( N/=nfcc_com ) THEN
              kp = jp + nmnr
              jp = jp + 1
              pjp = ry*DPRVEC(M,A(1,nr),A(1,j))
              P(jp) = pjp
              P(kp) = -pjp
              kp = kp + 1
              P(kp) = ry*dot
              DO k = 1, m2
                l = m2 + k
                A(k,j) = A(k,j) - pjp*A(l,nr)
                A(l,j) = A(l,j) + pjp*A(k,nr)
              END DO
              P(jq) = P(jq) - pjp*(pjp/ry)
            END IF
            !
            !                          TEST FOR CANCELLATION IN RECURRENCE RELATION
            !
            IF( P(jq)<=S(j)*sru_com ) P(jq) = NORM2(A(1:M,j))**2
            IF( P(jq)>y ) THEN
              y = P(jq)
              ix = j
            END IF
          END DO
          IF( N/=nfcc_com ) jp = kp
          !                       ************************************************
          IF( indpvt_com==1 ) ix = ip1
          !
          !                       RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH
          !                       SCALAR PRODUCT
          !
          y = NORM2(A(1:M,ix))**2
          !           ............EXIT
          IF( y<=eps_com*S(ix) ) GOTO 50
          Wcnd = MIN(Wcnd,y/S(ix))
        END IF
        !
        !                    COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR
        !                    SOLUTION
        !
        !                 ...EXIT
        IF( Inhomo==1 ) THEN
          lr = nr
          IF( N/=nfcc_com ) lr = 2*nr - 1
          W(lr) = DOT_PRODUCT(A(1:M,nr),V(1:M))*ry
          DO i = 1, M
            V(i) = V(i) - W(lr)*A(i,nr)
          END DO
          !                 ...EXIT
          IF( N/=nfcc_com ) THEN
            lr = 2*nr
            W(lr) = ry*DPRVEC(M,V,A(1,nr))
            DO k = 1, m2
              l = m2 + k
              V(k) = V(k) + W(lr)*A(l,nr)
              V(l) = V(l) - W(lr)*A(k,nr)
            END DO
          END IF
        END IF
      END DO
      !              *********************************************************
      !
      !                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
      !
      !        ......EXIT
      IF( Inhomo/=1 ) RETURN
      IF( (N>1) .AND. (S(np1)<1._SP) ) RETURN
      vnorm = NORM2(V(1:M))**2
      IF( S(np1)/=0._DP ) Wcnd = MIN(Wcnd,vnorm/S(np1))
      !        ......EXIT
      IF( vnorm>=eps_com*S(np1) ) RETURN
    END IF
    50  Iflag = 2
    Wcnd = eps_com
  ELSE
    Iflag = 1
  END IF
  !
  RETURN
END SUBROUTINE DMGSBV