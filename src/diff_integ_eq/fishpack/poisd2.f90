!** POISD2
SUBROUTINE POISD2(Mr,Nr,Istag,Ba,Bb,Bc,Q,Idimq,B,W,D,Tcos,P)
  !> Subsidiary to GENBUN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (POISD2-S, CMPOSD-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Subroutine to solve Poisson's equation for Dirichlet boundary
  !     conditions.
  !
  !     ISTAG = 1 if the last diagonal block is the matrix A.
  !     ISTAG = 2 if the last diagonal block is the matrix A+I.
  !
  !***
  ! **See also:**  GENBUN
  !***
  ! **Routines called:**  COSGEN, S1MERG, TRIX

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920130  Modified to use merge routine S1MERG rather than deleted
  !           routine MERGE.  (WRB)
  USE data_handling, ONLY : S1MERG
  INTEGER :: Idimq, Istag, Mr, Nr
  REAL(SP) :: B(Mr), Ba(Mr), Bb(Mr), Bc(Mr), D(Mr), P(:), Q(Idimq,Nr), Tcos(4*Nr), W(Mr)
  INTEGER :: i, ideg, ip, ip1, ipstor, irreg, j, jdeg, jm1, jm2, jm3, jp1, jp2, &
    jp3, jsh, jsp, jst, jstsav, kr, krpi, l, lr, m, n, nodd, noddpr, nun
  REAL(SP) :: fi, t
  !* FIRST EXECUTABLE STATEMENT  POISD2
  m = Mr
  n = Nr
  jsh = 0
  fi = 1._SP/Istag
  ip = -m
  ipstor = 0
  IF( Istag==2 ) THEN
    kr = 1
    jstsav = 1
    irreg = 2
    IF( n>1 ) GOTO 100
    Tcos(1) = -1._SP
  ELSE
    kr = 0
    irreg = 1
    IF( n>1 ) GOTO 100
    Tcos(1) = 0._SP
  END IF
  DO i = 1, m
    B(i) = Q(i,1)
  END DO
  CALL TRIX(1,0,m,Ba,Bb,Bc,B,Tcos,D,W)
  DO i = 1, m
    Q(i,1) = B(i)
  END DO
  !
  !     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
  !
  W(1) = ipstor
  RETURN
  100  lr = 0
  DO i = 1, m
    P(i) = 0._SP
  END DO
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
  IF( nodd==1 ) THEN
    jsp = jsp - jst
    IF( irreg/=1 ) jsp = jsp - l
  ELSE
    jsp = jsp - l
  END IF
  !
  !     REGULAR REDUCTION
  !
  CALL COSGEN(jst,1,0.5_SP,0._SP,Tcos)
  IF( l<=jsp ) THEN
    DO j = l, jsp, l
      jm1 = j - jsh
      jp1 = j + jsh
      jm2 = j - jst
      jp2 = j + jst
      jm3 = jm2 - jsh
      jp3 = jp2 + jsh
      IF( jst/=1 ) THEN
        DO i = 1, m
          t = Q(i,j) - Q(i,jm1) - Q(i,jp1) + Q(i,jm2) + Q(i,jp2)
          B(i) = t + Q(i,j) - Q(i,jm3) - Q(i,jp3)
          Q(i,j) = t
        END DO
      ELSE
        DO i = 1, m
          B(i) = 2._SP*Q(i,j)
          Q(i,j) = Q(i,jm2) + Q(i,jp2)
        END DO
      END IF
      CALL TRIX(jst,0,m,Ba,Bb,Bc,B,Tcos,D,W)
      DO i = 1, m
        Q(i,j) = Q(i,j) + B(i)
      END DO
    END DO
  END IF
  !
  !     REDUCTION FOR LAST UNKNOWN
  !
  IF( nodd==2 ) THEN
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
    IF( irreg==2 ) THEN
      CALL COSGEN(kr,jstsav,0._SP,fi,Tcos)
      CALL COSGEN(lr,jstsav,0._SP,fi,Tcos(kr+1))
      ideg = kr
      kr = kr + jst
    ELSE
      jstsav = jst
      ideg = jst
      kr = l
    END IF
    IF( jst/=1 ) THEN
      DO i = 1, m
        B(i) = Q(i,j) + 0.5_SP*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3))
      END DO
      IF( irreg/=2 ) THEN
        DO i = 1, m
          Q(i,j) = Q(i,jm2) + 0.5_SP*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
        END DO
        irreg = 2
      ELSEIF( noddpr==2 ) THEN
        DO i = 1, m
          Q(i,j) = Q(i,jm2) + Q(i,j) - Q(i,jm1)
        END DO
      ELSE
        DO i = 1, m
          ip1 = ip + i
          Q(i,j) = Q(i,jm2) + P(ip1)
        END DO
        ip = ip - m
      END IF
    ELSE
      irreg = 2
      DO i = 1, m
        B(i) = Q(i,j)
        Q(i,j) = Q(i,jm2)
      END DO
    END IF
    CALL TRIX(ideg,lr,m,Ba,Bb,Bc,B,Tcos,D,W)
    DO i = 1, m
      Q(i,j) = Q(i,j) + B(i)
    END DO
  ELSE
    IF( irreg==1 ) GOTO 300
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
    IF( Istag/=1 ) THEN
      IF( jst==1 ) THEN
        DO i = 1, m
          B(i) = Q(i,j)
          Q(i,j) = 0._SP
        END DO
        GOTO 250
      END IF
    END IF
    IF( noddpr==2 ) THEN
      DO i = 1, m
        B(i) = 0.5_SP*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3)) + Q(i,jp2) - Q(i,jp1) + Q(i,j)
      END DO
    ELSE
      DO i = 1, m
        ip1 = ip + i
        B(i) = 0.5_SP*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3)) + P(ip1) + Q(i,j)
      END DO
    END IF
    DO i = 1, m
      Q(i,j) = 0.5_SP*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
    END DO
    250  CALL TRIX(jst,0,m,Ba,Bb,Bc,B,Tcos,D,W)
    ip = ip + m
    ipstor = MAX(ipstor,ip+m)
    DO i = 1, m
      ip1 = ip + i
      P(ip1) = Q(i,j) + B(i)
      B(i) = Q(i,jp2) + P(ip1)
    END DO
    IF( lr/=0 ) THEN
      CALL COSGEN(lr,jstsav,0._SP,fi,Tcos(jst+1))
      CALL S1MERG(Tcos,0,jst,jst,lr,kr)
    ELSE
      DO i = 1, jst
        krpi = kr + i
        Tcos(krpi) = Tcos(i)
      END DO
    END IF
    CALL COSGEN(kr,jstsav,0._SP,fi,Tcos)
    CALL TRIX(kr,kr,m,Ba,Bb,Bc,B,Tcos,D,W)
    DO i = 1, m
      ip1 = ip + i
      Q(i,j) = Q(i,jm2) + B(i) + P(ip1)
    END DO
    lr = kr
    kr = kr + l
  END IF
  300  nun = nun/2
  noddpr = nodd
  jsh = jst
  jst = 2*jst
  IF( nun>=2 ) GOTO 200
  !
  !     START SOLUTION.
  !
  j = jsp
  DO i = 1, m
    B(i) = Q(i,j)
  END DO
  IF( irreg==2 ) THEN
    kr = lr + jst
    CALL COSGEN(kr,jstsav,0._SP,fi,Tcos)
    CALL COSGEN(lr,jstsav,0._SP,fi,Tcos(kr+1))
    ideg = kr
  ELSE
    CALL COSGEN(jst,1,0.5_SP,0._SP,Tcos)
    ideg = jst
  END IF
  CALL TRIX(ideg,lr,m,Ba,Bb,Bc,B,Tcos,D,W)
  jm1 = j - jsh
  jp1 = j + jsh
  IF( irreg/=2 ) THEN
    DO i = 1, m
      Q(i,j) = 0.5_SP*(Q(i,j)-Q(i,jm1)-Q(i,jp1)) + B(i)
    END DO
  ELSEIF( noddpr==2 ) THEN
    DO i = 1, m
      Q(i,j) = Q(i,j) - Q(i,jm1) + B(i)
    END DO
  ELSE
    DO i = 1, m
      ip1 = ip + i
      Q(i,j) = P(ip1) + B(i)
    END DO
    ip = ip - m
  END IF
  DO
    !
    !     START BACK SUBSTITUTION.
    !
    jst = jst/2
    jsh = jst/2
    nun = 2*nun
    IF( nun>n ) THEN
      W(1) = ipstor
      EXIT
    ELSE
      DO j = jst, n, l
        jm1 = j - jsh
        jp1 = j + jsh
        jm2 = j - jst
        jp2 = j + jst
        IF( j>jst ) THEN
          IF( jp2>n ) THEN
            DO i = 1, m
              B(i) = Q(i,j) + Q(i,jm2)
            END DO
            IF( jst<jstsav ) irreg = 1
            IF( irreg==1 ) GOTO 310
            IF( irreg==2 ) THEN
              IF( j+l>n ) lr = lr - jst
              kr = jst + lr
              CALL COSGEN(kr,jstsav,0._SP,fi,Tcos)
              CALL COSGEN(lr,jstsav,0._SP,fi,Tcos(kr+1))
              ideg = kr
              jdeg = lr
              GOTO 320
            END IF
          END IF
          DO i = 1, m
            B(i) = Q(i,j) + Q(i,jm2) + Q(i,jp2)
          END DO
        ELSE
          DO i = 1, m
            B(i) = Q(i,j) + Q(i,jp2)
          END DO
        END IF
        310  CALL COSGEN(jst,1,0.5_SP,0._SP,Tcos)
        ideg = jst
        jdeg = 0
        320  CALL TRIX(ideg,jdeg,m,Ba,Bb,Bc,B,Tcos,D,W)
        IF( jst>1 ) THEN
          IF( jp2>n ) THEN
            IF( irreg/=1 ) THEN
              IF( j+jsh>n ) THEN
                DO i = 1, m
                  Q(i,j) = B(i) + Q(i,j) - Q(i,jm1)
                END DO
              ELSE
                DO i = 1, m
                  ip1 = ip + i
                  Q(i,j) = B(i) + P(ip1)
                END DO
                ip = ip - m
              END IF
              CYCLE
            END IF
          END IF
          DO i = 1, m
            Q(i,j) = 0.5_SP*(Q(i,j)-Q(i,jm1)-Q(i,jp1)) + B(i)
          END DO
        ELSE
          DO i = 1, m
            Q(i,j) = B(i)
          END DO
        END IF
      END DO
      l = l/2
    END IF
  END DO
  RETURN
END SUBROUTINE POISD2
