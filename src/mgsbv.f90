!DECK MGSBV
SUBROUTINE MGSBV(M,N,A,Ia,Niv,Iflag,S,P,Ip,Inhomo,V,W,Wcnd)
  IMPLICIT NONE
  REAL A, AE, dot, EPS, FOUru, P, pjp, PRVEC, psave, RE, ry, S, &
    SDOT, SQOvfl, SRU, sv, t, TOL, TWOu, URO
  REAL V, vl, vnorm, W, Wcnd, y
  INTEGER i, Ia, ICOco, Iflag, INDpvt, Inhomo, INTeg, Ip, ip1, ix, &
    iz, j, jk, jp, jq, jy, jz, k, kd, kj
  INTEGER kp, l, lix, LPAr, lr, M, m2, MXNon, N, NDIsk, NEQ, &
    NEQivp, NFCc, NIC, Niv, nivn, nmnr, nn, NOPg, np1
  INTEGER NPS, nr, nrm1, NTApe, NTP, NUMort, NXPts
  !***BEGIN PROLOGUE  MGSBV
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (MGSBV-S, DMGSBV-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  ! Orthogonalize a set of N real vectors and determine their rank
  !
  ! **********************************************************************
  ! INPUT
  ! **********************************************************************
  !   M = Dimension of vectors
  !   N = No. of vectors
  !   A = Array whose first N cols contain the vectors
  !   IA = First dimension of array A (col length)
  !   NIV = Number of independent vectors needed
  !   INHOMO = 1 Corresponds to having a non-zero particular solution
  !   V = Particular solution vector (not included in the pivoting)
  !   INDPVT = 1 Means pivoting will not be used
  !
  ! **********************************************************************
  ! OUTPUT
  ! **********************************************************************
  !   NIV = No. of linear independent vectors in input set
  !     A = Matrix whose first NIV cols. contain NIV orthogonal vectors
  !         which span the vector space determined by the input vectors
  !   IFLAG
  !          = 0 success
  !          = 1 incorrect input
  !          = 2 rank of new vectors less than N
  !   P = Decomposition matrix.  P is upper triangular and
  !             (old vectors) = (new vectors) * P.
  !         The old vectors will be reordered due to pivoting
  !         The dimension of p must be .GE. N*(N+1)/2.
  !             (  N*(2*N+1) when N .NE. NFCC )
  !   IP = Pivoting vector. The dimension of IP must be .GE. N.
  !             (  2*N when N .NE. NFCC )
  !   S = Square of norms of incoming vectors
  !   V = Vector which is orthogonal to the vectors of A
  !   W = Orthogonalization information for the vector V
  !   WCND = Worst case (smallest) norm decrement value of the
  !          vectors being orthogonalized  (represents a test
  !          for linear dependence of the vectors)
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  PRVEC, SDOT
  !***COMMON BLOCKS    ML18JR, ML5MCO
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  MGSBV
  !
  DIMENSION A(Ia,*), V(*), W(*), P(*), Ip(*), S(*)
  !
  !
  COMMON /ML18JR/ AE, RE, TOL, NXPts, NIC, NOPg, MXNon, NDIsk, &
    NTApe, NEQ, INDpvt, INTeg, NPS, NTP, NEQivp, &
    NUMort, NFCc, ICOco
  !
  COMMON /ML5MCO/ URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  !
  !***FIRST EXECUTABLE STATEMENT  MGSBV
  IF ( M>0.AND.N>0.AND.Ia>=M ) THEN
    !
    jp = 0
    Iflag = 0
    np1 = N + 1
    y = 0.0
    m2 = M/2
    !
    !     CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH FOR
    !     VECTOR WITH LARGEST MAGNITUDE
    !
    j = 0
    DO i = 1, N
      vl = SDOT(M,A(1,i),1,A(1,i),1)
      S(i) = vl
      IF ( N/=NFCc ) THEN
        j = 2*i - 1
        P(j) = vl
        Ip(j) = j
      ENDIF
      j = j + 1
      P(j) = vl
      Ip(j) = j
      IF ( vl>y ) THEN
        y = vl
        ix = i
      ENDIF
    ENDDO
    IF ( INDpvt==1 ) THEN
      ix = 1
      y = P(1)
    ENDIF
    lix = ix
    IF ( N/=NFCc ) lix = 2*ix - 1
    P(lix) = P(1)
    S(np1) = 0.
    IF ( Inhomo==1 ) S(np1) = SDOT(M,V,1,V,1)
    Wcnd = 1.
    nivn = Niv
    Niv = 0
    !
    IF ( y/=0.0 ) THEN
      ! **********************************************************************
      DO nr = 1, N
        IF ( nivn==Niv ) EXIT
        Niv = nr
        IF ( ix/=nr ) THEN
          !
          !     PIVOTING OF COLUMNS OF P MATRIX
          !
          nn = N
          lix = ix
          lr = nr
          IF ( N/=NFCc ) THEN
            nn = NFCc
            lix = 2*ix - 1
            lr = 2*nr - 1
          ENDIF
          IF ( nr/=1 ) THEN
            kd = lix - lr
            kj = lr
            nrm1 = lr - 1
            DO j = 1, nrm1
              psave = P(kj)
              jk = kj + kd
              P(kj) = P(jk)
              P(jk) = psave
              kj = kj + nn - j
            ENDDO
            jy = jk + nmnr
            jz = jy - kd
            P(jy) = P(jz)
          ENDIF
          iz = Ip(lix)
          Ip(lix) = Ip(lr)
          Ip(lr) = iz
          sv = S(ix)
          S(ix) = S(nr)
          S(nr) = sv
          IF ( N/=NFCc ) THEN
            IF ( nr/=1 ) THEN
              kj = lr + 1
              DO k = 1, nrm1
                psave = P(kj)
                jk = kj + kd
                P(kj) = P(jk)
                P(jk) = psave
                kj = kj + NFCc - k
              ENDDO
            ENDIF
            iz = Ip(lix+1)
            Ip(lix+1) = Ip(lr+1)
            Ip(lr+1) = iz
          ENDIF
          !
          !     PIVOTING OF COLUMNS OF VECTORS
          !
          DO l = 1, M
            t = A(l,ix)
            A(l,ix) = A(l,nr)
            A(l,nr) = t
          ENDDO
        ENDIF
        !
        !     CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL VECTOR
        !
        jp = jp + 1
        P(jp) = y
        ry = 1.0/y
        nmnr = N - nr
        IF ( N/=NFCc ) THEN
          nmnr = NFCc - (2*nr-1)
          jp = jp + 1
          P(jp) = 0.
          kp = jp + nmnr
          P(kp) = y
        ENDIF
        IF ( nr/=N.AND.nivn/=Niv ) THEN
          !
          !    CALCULATE ORTHOGONAL PROJECTION VECTORS AND SEARCH FOR LARGEST NORM
          !
          y = 0.0
          ip1 = nr + 1
          ix = ip1
          !     ****************************************
          DO j = ip1, N
            dot = SDOT(M,A(1,nr),1,A(1,j),1)
            jp = jp + 1
            jq = jp + nmnr
            IF ( N/=NFCc ) jq = jq + nmnr - 1
            P(jq) = P(jp) - dot*(dot*ry)
            P(jp) = dot*ry
            DO i = 1, M
              A(i,j) = A(i,j) - P(jp)*A(i,nr)
            ENDDO
            IF ( N/=NFCc ) THEN
              kp = jp + nmnr
              jp = jp + 1
              pjp = ry*PRVEC(M,A(1,nr),A(1,j))
              P(jp) = pjp
              P(kp) = -pjp
              kp = kp + 1
              P(kp) = ry*dot
              DO k = 1, m2
                l = m2 + k
                A(k,j) = A(k,j) - pjp*A(l,nr)
                A(l,j) = A(l,j) + pjp*A(k,nr)
              ENDDO
              P(jq) = P(jq) - pjp*(pjp/ry)
            ENDIF
            !
            !     TEST FOR CANCELLATION IN RECURRENCE RELATION
            !
            IF ( P(jq)<=S(j)*SRU ) P(jq) = SDOT(M,A(1,j),1,A(1,j),1)
            IF ( P(jq)>y ) THEN
              y = P(jq)
              ix = j
            ENDIF
          ENDDO
          IF ( N/=NFCc ) jp = kp
          !     ****************************************
          IF ( INDpvt==1 ) ix = ip1
          !
          !     RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH SCALAR PRODUCT
          !
          y = SDOT(M,A(1,ix),1,A(1,ix),1)
          IF ( y<=EPS*S(ix) ) GOTO 50
          Wcnd = MIN(Wcnd,y/S(ix))
        ENDIF
        !
        !     COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR SOLUTION
        !
        IF ( Inhomo==1 ) THEN
          lr = nr
          IF ( N/=NFCc ) lr = 2*nr - 1
          W(lr) = SDOT(M,A(1,nr),1,V,1)*ry
          DO i = 1, M
            V(i) = V(i) - W(lr)*A(i,nr)
          ENDDO
          IF ( N/=NFCc ) THEN
            lr = 2*nr
            W(lr) = ry*PRVEC(M,V,A(1,nr))
            DO k = 1, m2
              l = m2 + k
              V(k) = V(k) + W(lr)*A(l,nr)
              V(l) = V(l) - W(lr)*A(k,nr)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      ! **********************************************************************
      !
      !     TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
      !
      IF ( Inhomo/=1 ) RETURN
      IF ( (N>1).AND.(S(np1)<1.0) ) RETURN
      vnorm = SDOT(M,V,1,V,1)
      IF ( S(np1)/=0. ) Wcnd = MIN(Wcnd,vnorm/S(np1))
      IF ( vnorm>=EPS*S(np1) ) RETURN
    ENDIF
    50     Iflag = 2
    Wcnd = EPS
  ELSE
    Iflag = 1
    RETURN
  ENDIF
END SUBROUTINE MGSBV
