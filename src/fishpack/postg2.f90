!DECK POSTG2
SUBROUTINE POSTG2(Nperod,N,M,A,Bb,C,Idimq,Q,B,B2,B3,W,W2,W3,D,Tcos,P)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  POSTG2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to POISTG
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (POSTG2-S)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Subroutine to solve Poisson's equation on a staggered grid.
  !
  !***SEE ALSO  POISTG
  !***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920130  Modified to use merge routine S1MERG rather than deleted
  !           routine MERGE.  (WRB)
  !***END PROLOGUE  POSTG2
  REAL A, B, B2, B3, Bb, C, D, fi, fnum, fnum2, P, Q, t, Tcos, W, W2, W3
  INTEGER i, i2r, i2rby2, Idimq, ii, ijump, ip, ipstor, j, jm1, &
    jm2, jm3, jp1, jp2, jp3, jr, jstart, jstep, jstop, k
  INTEGER k1, k2, k3, k4, kr, lr, M, mr, N, nlast, nlastp, np, &
    Nperod, nr, nrod, nrodpr
  DIMENSION A(*), Bb(*), C(*), Q(Idimq,*), B(*), B2(*), B3(*), W(*), &
    W2(*), W3(*), D(*), Tcos(*), k(4), P(*)
  EQUIVALENCE (k(1),k1)
  EQUIVALENCE (k(2),k2)
  EQUIVALENCE (k(3),k3)
  EQUIVALENCE (k(4),k4)
  !***FIRST EXECUTABLE STATEMENT  POSTG2
  np = Nperod
  fnum = 0.5*(np/3)
  fnum2 = 0.5*(np/2)
  mr = M
  ip = -mr
  ipstor = 0
  i2r = 1
  jr = 2
  nr = N
  nlast = N
  kr = 1
  lr = 0
  IF ( nr<=3 ) GOTO 200
  100  jr = 2*i2r
  nrod = 1
  IF ( (nr/2)*2==nr ) nrod = 0
  jstart = 1
  jstop = nlast - jr
  IF ( nrod==0 ) jstop = jstop - i2r
  i2rby2 = i2r/2
  IF ( jstop>=jstart ) THEN
    !
    !     REGULAR REDUCTION.
    !
    ijump = 1
    DO j = jstart, jstop, jr
      jp1 = j + i2rby2
      jp2 = j + i2r
      jp3 = jp2 + i2rby2
      jm1 = j - i2rby2
      jm2 = j - i2r
      jm3 = jm2 - i2rby2
      IF ( j/=1 ) THEN
        IF ( ijump/=2 ) THEN
          ijump = 2
          CALL COSGEN(i2r,1,0.5,0.0,Tcos)
        ENDIF
        IF ( i2r/=1 ) THEN
          DO i = 1, mr
            fi = Q(i,j)
            Q(i,j) = Q(i,j) - Q(i,jm1) - Q(i,jp1) + Q(i,jm2) + Q(i,jp2)
            B(i) = fi + Q(i,j) - Q(i,jm3) - Q(i,jp3)
          ENDDO
        ELSE
          DO i = 1, mr
            B(i) = 2.*Q(i,j)
            Q(i,j) = Q(i,jm2) + Q(i,jp2)
          ENDDO
        ENDIF
      ELSE
        CALL COSGEN(i2r,1,fnum,0.5,Tcos)
        IF ( i2r/=1 ) THEN
          DO i = 1, mr
            B(i) = Q(i,1) + 0.5*(Q(i,jp2)-Q(i,jp1)-Q(i,jp3))
            Q(i,1) = Q(i,jp2) + Q(i,1) - Q(i,jp1)
          ENDDO
        ELSE
          DO i = 1, mr
            B(i) = Q(i,1)
            Q(i,1) = Q(i,2)
          ENDDO
        ENDIF
      ENDIF
      CALL TRIX(i2r,0,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,j) = Q(i,j) + B(i)
      ENDDO
      !
      !     END OF REDUCTION FOR REGULAR UNKNOWNS.
      !
    ENDDO
    !
    !     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
    !
    j = jstop + jr
  ELSE
    j = jr
  ENDIF
  nlast = j
  jm1 = j - i2rby2
  jm2 = j - i2r
  jm3 = jm2 - i2rby2
  IF ( nrod==0 ) THEN
    !
    !     EVEN NUMBER OF UNKNOWNS
    !
    jp1 = j + i2rby2
    jp2 = j + i2r
    IF ( i2r/=1 ) THEN
      DO i = 1, mr
        B(i) = Q(i,j) + .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3))
      ENDDO
      IF ( nrodpr/=0 ) THEN
        DO i = 1, mr
          B(i) = B(i) + Q(i,jp2) - Q(i,jp1)
        ENDDO
      ELSE
        DO i = 1, mr
          ii = ip + i
          B(i) = B(i) + P(ii)
        ENDDO
      ENDIF
      CALL COSGEN(i2r,1,0.5,0.0,Tcos)
      CALL TRIX(i2r,0,mr,A,Bb,C,B,Tcos,D,W)
      ip = ip + mr
      ipstor = MAX(ipstor,ip+mr)
      DO i = 1, mr
        ii = ip + i
        P(ii) = B(i) + .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
        B(i) = P(ii) + Q(i,jp2)
      ENDDO
      IF ( lr==0 ) THEN
        DO i = 1, i2r
          ii = kr + i
          Tcos(ii) = Tcos(i)
        ENDDO
      ELSE
        CALL COSGEN(lr,1,fnum2,0.5,Tcos(i2r+1))
        CALL S1MERG(Tcos,0,i2r,i2r,lr,kr)
      ENDIF
      CALL COSGEN(kr,1,fnum2,0.5,Tcos)
      CALL TRIX(kr,kr,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        ii = ip + i
        Q(i,j) = Q(i,jm2) + P(ii) + B(i)
      ENDDO
    ELSE
      DO i = 1, mr
        B(i) = Q(i,j)
      ENDDO
      Tcos(1) = 0.
      CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
      ip = 0
      ipstor = mr
      DO i = 1, mr
        P(i) = B(i)
        B(i) = B(i) + Q(i,N)
      ENDDO
      Tcos(1) = -1. + 2*(np/2)
      Tcos(2) = 0.
      CALL TRIX(1,1,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,j) = Q(i,jm2) + P(i) + B(i)
      ENDDO
    ENDIF
    lr = kr
    kr = kr + jr
  ELSE
    !
    !     ODD NUMBER OF UNKNOWNS
    !
    IF ( i2r/=1 ) THEN
      DO i = 1, mr
        B(i) = Q(i,j) + .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3))
      ENDDO
      IF ( nrodpr/=0 ) THEN
        DO i = 1, mr
          Q(i,j) = Q(i,j) - Q(i,jm1) + Q(i,jm2)
        ENDDO
      ELSE
        DO i = 1, mr
          ii = ip + i
          Q(i,j) = Q(i,jm2) + P(ii)
        ENDDO
        ip = ip - mr
      ENDIF
      IF ( lr/=0 ) CALL COSGEN(lr,1,fnum2,0.5,Tcos(kr+1))
    ELSE
      DO i = 1, mr
        B(i) = Q(i,j)
        Q(i,j) = Q(i,jm2)
      ENDDO
    ENDIF
    CALL COSGEN(kr,1,fnum2,0.5,Tcos)
    CALL TRIX(kr,lr,mr,A,Bb,C,B,Tcos,D,W)
    DO i = 1, mr
      Q(i,j) = Q(i,j) + B(i)
    ENDDO
    kr = kr + i2r
  ENDIF
  nr = (nlast-1)/jr + 1
  IF ( nr>3 ) THEN
    i2r = jr
    nrodpr = nrod
    GOTO 100
  ENDIF
  !
  !      BEGIN SOLUTION
  !
  200  j = 1 + jr
  jm1 = j - i2r
  jp1 = j + i2r
  jm2 = nlast - i2r
  IF ( nr==2 ) THEN
    !
    !     CASE OF GENERAL N AND NR = 2 .
    !
    DO i = 1, mr
      ii = ip + i
      B3(i) = 0.
      B(i) = Q(i,1) + P(ii)
      Q(i,1) = Q(i,1) - Q(i,jm1)
      B2(i) = Q(i,1) + Q(i,nlast)
    ENDDO
    k1 = kr + jr
    k2 = k1 + jr
    CALL COSGEN(jr-1,1,0.0,1.0,Tcos(k1+1))
    SELECT CASE (np)
      CASE (2)
        CALL COSGEN(kr+1,1,0.5,0.0,Tcos(k2))
      CASE DEFAULT
        Tcos(k2) = 2*np - 4
        CALL COSGEN(kr,1,0.0,1.0,Tcos(k2+1))
    END SELECT
    k4 = 1 - np/3
    CALL S1MERG(Tcos,k1,jr-k4,k2-k4,kr+k4,0)
    IF ( np==3 ) k1 = k1 - 1
    k2 = kr
    CALL COSGEN(kr,1,fnum2,0.5,Tcos(k1+1))
    k4 = k1 + kr
    CALL COSGEN(lr,1,fnum2,0.5,Tcos(k4+1))
    k3 = lr
    k4 = 0
    CALL TRI3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
    DO i = 1, mr
      B(i) = B(i) + B2(i)
    ENDDO
    IF ( np==3 ) THEN
      Tcos(1) = 2.
      CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
    ENDIF
    DO i = 1, mr
      Q(i,1) = Q(i,1) + B(i)
    ENDDO
  ELSEIF ( lr/=0 ) THEN
    !
    !     CASE OF GENERAL N WITH NR = 3 .
    !
    DO i = 1, mr
      B(i) = Q(i,1) - Q(i,jm1) + Q(i,j)
    ENDDO
    IF ( nrod/=0 ) THEN
      DO i = 1, mr
        B(i) = B(i) + Q(i,nlast) - Q(i,jm2)
      ENDDO
    ELSE
      DO i = 1, mr
        ii = ip + i
        B(i) = B(i) + P(ii)
      ENDDO
    ENDIF
    DO i = 1, mr
      t = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
      Q(i,j) = t
      B2(i) = Q(i,nlast) + t
      B3(i) = Q(i,1) + t
    ENDDO
    k1 = kr + 2*jr
    CALL COSGEN(jr-1,1,0.0,1.0,Tcos(k1+1))
    k2 = k1 + jr
    Tcos(k2) = 2*np - 4
    k4 = (np-1)*(3-np)
    k3 = k2 + 1 - k4
    CALL COSGEN(kr+jr+k4,1,k4/2.,1.-k4,Tcos(k3))
    k4 = 1 - np/3
    CALL S1MERG(Tcos,k1,jr-k4,k2-k4,kr+jr+k4,0)
    IF ( np==3 ) k1 = k1 - 1
    k2 = kr + jr
    k4 = k1 + k2
    CALL COSGEN(kr,1,fnum2,0.5,Tcos(k4+1))
    k3 = k4 + kr
    CALL COSGEN(jr,1,fnum,0.5,Tcos(k3+1))
    CALL S1MERG(Tcos,k4,kr,k3,jr,k1)
    k4 = k3 + jr
    CALL COSGEN(lr,1,fnum2,0.5,Tcos(k4+1))
    CALL S1MERG(Tcos,k3,jr,k4,lr,k1+k2)
    CALL COSGEN(kr,1,fnum2,0.5,Tcos(k3+1))
    k3 = kr
    k4 = kr
    CALL TRI3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
    DO i = 1, mr
      B(i) = B(i) + B2(i) + B3(i)
    ENDDO
    IF ( np==3 ) THEN
      Tcos(1) = 2.
      CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
    ENDIF
    DO i = 1, mr
      Q(i,j) = Q(i,j) + B(i)
      B(i) = Q(i,1) + Q(i,j)
    ENDDO
    CALL COSGEN(jr,1,fnum,0.5,Tcos)
    CALL TRIX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
    IF ( jr/=1 ) THEN
      DO i = 1, mr
        Q(i,1) = Q(i,1) - Q(i,jm1) + B(i)
      ENDDO
    ELSE
      DO i = 1, mr
        Q(i,1) = B(i)
      ENDDO
    ENDIF
  ELSEIF ( N/=3 ) THEN
    !
    !     CASE N = 2**P+1
    !
    DO i = 1, mr
      B(i) = Q(i,j) + Q(i,1) - Q(i,jm1) + Q(i,nlast) - Q(i,jm2)
    ENDDO
    SELECT CASE (np)
      CASE (2)
        DO i = 1, mr
          fi = (Q(i,j)-Q(i,jm1)-Q(i,jp1))/2.
          B2(i) = Q(i,1) + fi
          B3(i) = Q(i,nlast) + fi
        ENDDO
        k1 = nlast + jr - 1
        k2 = k1 + jr - 1
        CALL COSGEN(jr-1,1,0.0,1.0,Tcos(k1+1))
        CALL COSGEN(nlast,1,0.5,0.0,Tcos(k2+1))
        CALL S1MERG(Tcos,k1,jr-1,k2,nlast,0)
        k3 = k1 + nlast - 1
        k4 = k3 + jr
        CALL COSGEN(jr,1,0.5,0.5,Tcos(k3+1))
        CALL COSGEN(jr,1,0.0,0.5,Tcos(k4+1))
        CALL S1MERG(Tcos,k3,jr,k4,jr,k1)
        k2 = nlast - 1
        k3 = jr
        k4 = jr
      CASE DEFAULT
        DO i = 1, mr
          B2(i) = Q(i,1) + Q(i,nlast) + Q(i,j) - Q(i,jm1) - Q(i,jp1)
          B3(i) = 0.
        ENDDO
        k1 = nlast - 1
        k2 = nlast + jr - 1
        CALL COSGEN(jr-1,1,0.0,1.0,Tcos(nlast))
        Tcos(k2) = 2*np - 4
        CALL COSGEN(jr,1,0.5-fnum,0.5,Tcos(k2+1))
        k3 = (3-np)/2
        CALL S1MERG(Tcos,k1,jr-k3,k2-k3,jr+k3,0)
        k1 = k1 - 1 + k3
        CALL COSGEN(jr,1,fnum,0.5,Tcos(k1+1))
        k2 = jr
        k3 = 0
        k4 = 0
    END SELECT
    CALL TRI3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
    DO i = 1, mr
      B(i) = B(i) + B2(i) + B3(i)
    ENDDO
    IF ( np==3 ) THEN
      Tcos(1) = 2.
      CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
    ENDIF
    DO i = 1, mr
      Q(i,j) = B(i) + .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
      B(i) = Q(i,j) + Q(i,1)
    ENDDO
    CALL COSGEN(jr,1,fnum,0.5,Tcos)
    CALL TRIX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
    DO i = 1, mr
      Q(i,1) = Q(i,1) - Q(i,jm1) + B(i)
    ENDDO
  ELSE
    !
    !     CASE N = 3.
    !
    SELECT CASE (np)
      CASE (2)
        DO i = 1, mr
          B(i) = Q(i,2)
          B2(i) = Q(i,3)
          B3(i) = Q(i,1)
        ENDDO
        CALL COSGEN(3,1,0.5,0.0,Tcos)
        Tcos(4) = -1.
        Tcos(5) = 1.
        Tcos(6) = -1.
        Tcos(7) = 1.
        k1 = 3
        k2 = 2
        k3 = 1
        k4 = 1
      CASE DEFAULT
        DO i = 1, mr
          B(i) = Q(i,2)
          B2(i) = Q(i,1) + Q(i,3)
          B3(i) = 0.
        ENDDO
        SELECT CASE (np)
          CASE (1,2)
            Tcos(1) = -2.
            Tcos(2) = 1.
            Tcos(3) = -1.
            k1 = 2
          CASE DEFAULT
            Tcos(1) = -1.
            Tcos(2) = 1.
            k1 = 1
        END SELECT
        k2 = 1
        k3 = 0
        k4 = 0
    END SELECT
    CALL TRI3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
    DO i = 1, mr
      B(i) = B(i) + B2(i) + B3(i)
    ENDDO
    SELECT CASE (np)
      CASE (1,2)
      CASE DEFAULT
        Tcos(1) = 2.
        CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
    END SELECT
    DO i = 1, mr
      Q(i,2) = B(i)
      B(i) = Q(i,1) + B(i)
    ENDDO
    Tcos(1) = -1. + 4.*fnum
    CALL TRIX(1,0,mr,A,Bb,C,B,Tcos,D,W)
    DO i = 1, mr
      Q(i,1) = B(i)
    ENDDO
    jr = 1
    i2r = 0
  ENDIF
  !
  !     START BACK SUBSTITUTION.
  !
  300  j = nlast - jr
  DO i = 1, mr
    B(i) = Q(i,nlast) + Q(i,j)
  ENDDO
  jm2 = nlast - i2r
  IF ( jr==1 ) THEN
    DO i = 1, mr
      Q(i,nlast) = 0.
    ENDDO
  ELSEIF ( nrod/=0 ) THEN
    DO i = 1, mr
      Q(i,nlast) = Q(i,nlast) - Q(i,jm2)
    ENDDO
  ELSE
    DO i = 1, mr
      ii = ip + i
      Q(i,nlast) = P(ii)
    ENDDO
    ip = ip - mr
  ENDIF
  CALL COSGEN(kr,1,fnum2,0.5,Tcos)
  CALL COSGEN(lr,1,fnum2,0.5,Tcos(kr+1))
  CALL TRIX(kr,lr,mr,A,Bb,C,B,Tcos,D,W)
  DO i = 1, mr
    Q(i,nlast) = Q(i,nlast) + B(i)
  ENDDO
  nlastp = nlast
  DO
    jstep = jr
    jr = i2r
    i2r = i2r/2
    IF ( jr==0 ) THEN
      !
      !     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
      !
      W(1) = ipstor
      EXIT
    ELSE
      jstart = 1 + jr
      kr = kr - jr
      IF ( nlast+jr>N ) THEN
        jstop = nlast - jr
      ELSE
        kr = kr - jr
        nlast = nlast + jr
        jstop = nlast - jstep
      ENDIF
      lr = kr - jr
      CALL COSGEN(jr,1,0.5,0.0,Tcos)
      DO j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr
        IF ( j/=jr ) THEN
          DO i = 1, mr
            B(i) = Q(i,j) + Q(i,jm2) + Q(i,jp2)
          ENDDO
        ELSE
          DO i = 1, mr
            B(i) = Q(i,j) + Q(i,jp2)
          ENDDO
        ENDIF
        IF ( jr/=1 ) THEN
          jm1 = j - i2r
          jp1 = j + i2r
          DO i = 1, mr
            Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
          ENDDO
        ELSE
          DO i = 1, mr
            Q(i,j) = 0.
          ENDDO
        ENDIF
        CALL TRIX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,j) = Q(i,j) + B(i)
        ENDDO
      ENDDO
      nrod = 1
      IF ( nlast+i2r<=N ) nrod = 0
      IF ( nlastp/=nlast ) GOTO 300
    ENDIF
  ENDDO
END SUBROUTINE POSTG2
