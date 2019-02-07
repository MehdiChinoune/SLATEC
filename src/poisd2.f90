!*==POISD2.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK POISD2
SUBROUTINE POISD2(Mr,Nr,Istag,Ba,Bb,Bc,Q,Idimq,B,W,D,Tcos,P)
  IMPLICIT NONE
  !*--POISD25
  !*** Start of declarations inserted by SPAG
  REAL B, Ba, Bb, Bc, D, fi, P, Q, t, Tcos, W
  INTEGER i, ideg, Idimq, ip, ip1, ipstor, irreg, Istag, j, jdeg, &
    jm1, jm2, jm3, jp1, jp2, jp3, jsh, jsp, jst, jstsav
  INTEGER kr, krpi, l, lr, m, Mr, n, nodd, noddpr, Nr, nun
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  POISD2
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to GENBUN
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (POISD2-S, CMPOSD-C)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     Subroutine to solve Poisson's equation for Dirichlet boundary
  !     conditions.
  !
  !     ISTAG = 1 if the last diagonal block is the matrix A.
  !     ISTAG = 2 if the last diagonal block is the matrix A+I.
  !
  !***SEE ALSO  GENBUN
  !***ROUTINES CALLED  COSGEN, S1MERG, TRIX
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920130  Modified to use merge routine S1MERG rather than deleted
  !           routine MERGE.  (WRB)
  !***END PROLOGUE  POISD2
  !
  DIMENSION Q(Idimq,*), Ba(*), Bb(*), Bc(*), Tcos(*), B(*), D(*), &
    W(*), P(*)
  !***FIRST EXECUTABLE STATEMENT  POISD2
  m = Mr
  n = Nr
  jsh = 0
  fi = 1./Istag
  ip = -m
  ipstor = 0
  IF ( Istag==2 ) THEN
    kr = 1
    jstsav = 1
    irreg = 2
    IF ( n>1 ) GOTO 100
    Tcos(1) = -1.
  ELSE
    kr = 0
    irreg = 1
    IF ( n>1 ) GOTO 100
    Tcos(1) = 0.
  ENDIF
  DO i = 1, m
    B(i) = Q(i,1)
  ENDDO
  CALL TRIX(1,0,m,Ba,Bb,Bc,B,Tcos,D,W)
  DO i = 1, m
    Q(i,1) = B(i)
  ENDDO
  !
  !     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
  !
  W(1) = ipstor
  GOTO 99999
  100  lr = 0
  DO i = 1, m
    P(i) = 0.
  ENDDO
  nun = n
  jst = 1
  jsp = n
  !
  !     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
  !
  200  l = 2*jst
  nodd = 2 - 2*((nun+1)/2) + nun
  !
  !     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
  !
  IF ( nodd==1 ) THEN
    jsp = jsp - jst
    IF ( irreg/=1 ) jsp = jsp - l
  ELSE
    jsp = jsp - l
  ENDIF
  !
  !     REGULAR REDUCTION
  !
  CALL COSGEN(jst,1,0.5,0.0,Tcos)
  IF ( l<=jsp ) THEN
    DO j = l, jsp, l
      jm1 = j - jsh
      jp1 = j + jsh
      jm2 = j - jst
      jp2 = j + jst
      jm3 = jm2 - jsh
      jp3 = jp2 + jsh
      IF ( jst/=1 ) THEN
        DO i = 1, m
          t = Q(i,j) - Q(i,jm1) - Q(i,jp1) + Q(i,jm2) + Q(i,jp2)
          B(i) = t + Q(i,j) - Q(i,jm3) - Q(i,jp3)
          Q(i,j) = t
        ENDDO
      ELSE
        DO i = 1, m
          B(i) = 2.*Q(i,j)
          Q(i,j) = Q(i,jm2) + Q(i,jp2)
        ENDDO
      ENDIF
      CALL TRIX(jst,0,m,Ba,Bb,Bc,B,Tcos,D,W)
      DO i = 1, m
        Q(i,j) = Q(i,j) + B(i)
      ENDDO
    ENDDO
  ENDIF
  !
  !     REDUCTION FOR LAST UNKNOWN
  !
  IF ( nodd==2 ) THEN
    !
    !     EVEN NUMBER OF UNKNOWNS
    !
    jsp = jsp + l
    j = jsp
    jm1 = j - jsh
    jp1 = j + jsh
    jm2 = j - jst
    jp2 = j + jst
    jm3 = jm2 - jsh
    IF ( irreg==2 ) THEN
      CALL COSGEN(kr,jstsav,0.0,fi,Tcos)
      CALL COSGEN(lr,jstsav,0.0,fi,Tcos(kr+1))
      ideg = kr
      kr = kr + jst
    ELSE
      jstsav = jst
      ideg = jst
      kr = l
    ENDIF
    IF ( jst/=1 ) THEN
      DO i = 1, m
        B(i) = Q(i,j) + .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3))
      ENDDO
      IF ( irreg/=2 ) THEN
        DO i = 1, m
          Q(i,j) = Q(i,jm2) + .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
        ENDDO
        irreg = 2
      ELSEIF ( noddpr==2 ) THEN
        DO i = 1, m
          Q(i,j) = Q(i,jm2) + Q(i,j) - Q(i,jm1)
        ENDDO
      ELSE
        DO i = 1, m
          ip1 = ip + i
          Q(i,j) = Q(i,jm2) + P(ip1)
        ENDDO
        ip = ip - m
      ENDIF
    ELSE
      irreg = 2
      DO i = 1, m
        B(i) = Q(i,j)
        Q(i,j) = Q(i,jm2)
      ENDDO
    ENDIF
    CALL TRIX(ideg,lr,m,Ba,Bb,Bc,B,Tcos,D,W)
    DO i = 1, m
      Q(i,j) = Q(i,j) + B(i)
    ENDDO
  ELSE
    IF ( irreg==1 ) GOTO 300
    !
    !     ODD NUMBER OF UNKNOWNS
    !
    jsp = jsp + l
    j = jsp
    jm1 = j - jsh
    jp1 = j + jsh
    jm2 = j - jst
    jp2 = j + jst
    jm3 = jm2 - jsh
    IF ( Istag/=1 ) THEN
      IF ( jst==1 ) THEN
        DO i = 1, m
          B(i) = Q(i,j)
          Q(i,j) = 0.
        ENDDO
        GOTO 250
      ENDIF
    ENDIF
    IF ( noddpr==2 ) THEN
      DO i = 1, m
        B(i) = .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3)) + Q(i,jp2) - Q(i,jp1)&
          + Q(i,j)
      ENDDO
    ELSE
      DO i = 1, m
        ip1 = ip + i
        B(i) = .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3)) + P(ip1) + Q(i,j)
      ENDDO
    ENDIF
    DO i = 1, m
      Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
    ENDDO
    250    CALL TRIX(jst,0,m,Ba,Bb,Bc,B,Tcos,D,W)
    ip = ip + m
    ipstor = MAX(ipstor,ip+m)
    DO i = 1, m
      ip1 = ip + i
      P(ip1) = Q(i,j) + B(i)
      B(i) = Q(i,jp2) + P(ip1)
    ENDDO
    IF ( lr/=0 ) THEN
      CALL COSGEN(lr,jstsav,0.,fi,Tcos(jst+1))
      CALL S1MERG(Tcos,0,jst,jst,lr,kr)
    ELSE
      DO i = 1, jst
        krpi = kr + i
        Tcos(krpi) = Tcos(i)
      ENDDO
    ENDIF
    CALL COSGEN(kr,jstsav,0.0,fi,Tcos)
    CALL TRIX(kr,kr,m,Ba,Bb,Bc,B,Tcos,D,W)
    DO i = 1, m
      ip1 = ip + i
      Q(i,j) = Q(i,jm2) + B(i) + P(ip1)
    ENDDO
    lr = kr
    kr = kr + l
  ENDIF
  300  nun = nun/2
  noddpr = nodd
  jsh = jst
  jst = 2*jst
  IF ( nun>=2 ) GOTO 200
  !
  !     START SOLUTION.
  !
  j = jsp
  DO i = 1, m
    B(i) = Q(i,j)
  ENDDO
  IF ( irreg==2 ) THEN
    kr = lr + jst
    CALL COSGEN(kr,jstsav,0.0,fi,Tcos)
    CALL COSGEN(lr,jstsav,0.0,fi,Tcos(kr+1))
    ideg = kr
  ELSE
    CALL COSGEN(jst,1,0.5,0.0,Tcos)
    ideg = jst
  ENDIF
  CALL TRIX(ideg,lr,m,Ba,Bb,Bc,B,Tcos,D,W)
  jm1 = j - jsh
  jp1 = j + jsh
  IF ( irreg/=2 ) THEN
    DO i = 1, m
      Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1)) + B(i)
    ENDDO
  ELSEIF ( noddpr==2 ) THEN
    DO i = 1, m
      Q(i,j) = Q(i,j) - Q(i,jm1) + B(i)
    ENDDO
  ELSE
    DO i = 1, m
      ip1 = ip + i
      Q(i,j) = P(ip1) + B(i)
    ENDDO
    ip = ip - m
  ENDIF
  DO
    !
    !     START BACK SUBSTITUTION.
    !
    jst = jst/2
    jsh = jst/2
    nun = 2*nun
    IF ( nun>n ) THEN
      W(1) = ipstor
      EXIT
    ELSE
      DO j = jst, n, l
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        IF ( j>jst ) THEN
          IF ( jp2>n ) THEN
            DO i = 1, m
              B(i) = Q(i,j) + Q(i,jm2)
            ENDDO
            IF ( jst<jstsav ) irreg = 1
            IF ( irreg==1 ) GOTO 310
            IF ( irreg==2 ) THEN
              IF ( j+l>n ) lr = lr - jst
              kr = jst + lr
              CALL COSGEN(kr,jstsav,0.0,fi,Tcos)
              CALL COSGEN(lr,jstsav,0.0,fi,Tcos(kr+1))
              ideg = kr
              jdeg = lr
              GOTO 320
            ENDIF
          ENDIF
          DO i = 1, m
            B(i) = Q(i,j) + Q(i,jm2) + Q(i,jp2)
          ENDDO
        ELSE
          DO i = 1, m
            B(i) = Q(i,j) + Q(i,jp2)
          ENDDO
        ENDIF
        310        CALL COSGEN(jst,1,0.5,0.0,Tcos)
        ideg = jst
        jdeg = 0
        320        CALL TRIX(ideg,jdeg,m,Ba,Bb,Bc,B,Tcos,D,W)
        IF ( jst>1 ) THEN
          IF ( jp2>n ) THEN
            IF ( irreg/=1 ) THEN
              IF ( j+jsh>n ) THEN
                DO i = 1, m
                  Q(i,j) = B(i) + Q(i,j) - Q(i,jm1)
                ENDDO
              ELSE
                DO i = 1, m
                  ip1 = ip + i
                  Q(i,j) = B(i) + P(ip1)
                ENDDO
                ip = ip - m
              ENDIF
              CYCLE
            ENDIF
          ENDIF
          DO i = 1, m
            Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1)) + B(i)
          ENDDO
        ELSE
          DO i = 1, m
            Q(i,j) = B(i)
          ENDDO
        ENDIF
      ENDDO
      l = l/2
    ENDIF
  ENDDO
  99999 CONTINUE
  END SUBROUTINE POISD2
