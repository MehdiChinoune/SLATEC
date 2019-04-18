!** CMPOSN
SUBROUTINE CMPOSN(M,N,Istag,Mixbnd,A,Bb,C,Q,Idimq,B,B2,B3,W,W2,W3,D,Tcos,P)
  !>
  !  Subsidiary to CMGNBN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      COMPLEX (POISN2-S, CMPOSN-C)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     Subroutine to solve Poisson's equation with Neumann boundary
  !     conditions.
  !
  !     ISTAG = 1 if the last diagonal block is A.
  !     ISTAG = 2 if the last diagonal block is A-I.
  !     MIXBND = 1 if have Neumann boundary conditions at both boundaries.
  !     MIXBND = 2 if have Neumann boundary conditions at bottom and
  !     Dirichlet condition at top.  (For this case, must have ISTAG = 1)
  !
  !***
  ! **See also:**  CMGNBN
  !***
  ! **Routines called:**  C1MERG, CMPCSG, CMPTR3, CMPTRX

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !   920130  Modified to use merge routine C1MERG rather than deleted
  !           routine CMPMRG.  (WRB)
  USE data_handling, ONLY : C1MERG
  REAL fden, fistag, fnum
  INTEGER i, i1, i2, i2r, i2rby2, Idimq, ii, ip, ipstor, Istag, j, jm1, jm2, jm3, &
    jp1, jp2, jp3, jr, jr2, jstart, jstep, jstop, k(4), kr, lr, M, Mixbnd, mr, N, &
    nlast, nlastp, nr, nrod, nrodpr
  COMPLEX A(*), Bb(*), C(*), Q(Idimq,*), B(*), B2(*), B3(*), W(*), W2(*), &
    W3(*), D(*), Tcos(*), P(*), fi, t
  !* FIRST EXECUTABLE STATEMENT  CMPOSN
  fistag = 3 - Istag
  fnum = 1./Istag
  fden = 0.5*(Istag-1)
  mr = M
  ip = -mr
  ipstor = 0
  i2r = 1
  jr = 2
  nr = N
  nlast = N
  kr = 1
  lr = 0
  IF ( Istag/=2 ) THEN
    DO i = 1, mr
      Q(i,N) = .5*Q(i,N)
    END DO
    IF ( Mixbnd==2 ) GOTO 100
  END IF
  IF ( N<=3 ) GOTO 200
  100  jr = 2*i2r
  nrod = 1
  IF ( (nr/2)*2==nr ) nrod = 0
  IF ( Mixbnd==2 ) THEN
    jstart = jr
    nrod = 1 - nrod
  ELSE
    jstart = 1
  END IF
  jstop = nlast - jr
  IF ( nrod==0 ) jstop = jstop - i2r
  CALL CMPCSG(i2r,1,0.5,0.0,Tcos)
  i2rby2 = i2r/2
  IF ( jstop>=jstart ) THEN
    !
    !     REGULAR REDUCTION.
    !
    DO j = jstart, jstop, jr
      jp1 = j + i2rby2
      jp2 = j + i2r
      jp3 = jp2 + i2rby2
      jm1 = j - i2rby2
      jm2 = j - i2r
      jm3 = jm2 - i2rby2
      IF ( j==1 ) THEN
        jm1 = jp1
        jm2 = jp2
        jm3 = jp3
      END IF
      IF ( i2r/=1 ) THEN
        DO i = 1, mr
          fi = Q(i,j)
          Q(i,j) = Q(i,j) - Q(i,jm1) - Q(i,jp1) + Q(i,jm2) + Q(i,jp2)
          B(i) = fi + Q(i,j) - Q(i,jm3) - Q(i,jp3)
        END DO
      ELSE
        IF ( j==1 ) jm2 = jp2
        DO i = 1, mr
          B(i) = 2.*Q(i,j)
          Q(i,j) = Q(i,jm2) + Q(i,jp2)
        END DO
      END IF
      CALL CMPTRX(i2r,0,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,j) = Q(i,j) + B(i)
      END DO
      !
      !     END OF REDUCTION FOR REGULAR UNKNOWNS.
      !
    END DO
    !
    !     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
    !
    j = jstop + jr
  ELSE
    j = jr
  END IF
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
      END DO
      IF ( nrodpr/=0 ) THEN
        DO i = 1, mr
          B(i) = B(i) + Q(i,jp2) - Q(i,jp1)
        END DO
      ELSE
        DO i = 1, mr
          ii = ip + i
          B(i) = B(i) + P(ii)
        END DO
      END IF
      CALL CMPTRX(i2r,0,mr,A,Bb,C,B,Tcos,D,W)
      ip = ip + mr
      ipstor = MAX(ipstor,ip+mr)
      DO i = 1, mr
        ii = ip + i
        P(ii) = B(i) + .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
        B(i) = P(ii) + Q(i,jp2)
      END DO
      IF ( lr==0 ) THEN
        DO i = 1, i2r
          ii = kr + i
          Tcos(ii) = Tcos(i)
        END DO
      ELSE
        CALL CMPCSG(lr,1,0.5,fden,Tcos(i2r+1))
        CALL C1MERG(Tcos,0,i2r,i2r,lr,kr)
      END IF
      CALL CMPCSG(kr,1,0.5,fden,Tcos)
      IF ( lr/=0 ) THEN
        CALL CMPTRX(kr,kr,mr,A,Bb,C,B,Tcos,D,W)
      ELSEIF ( Istag==1 ) THEN
        DO i = 1, mr
          B(i) = fistag*B(i)
        END DO
      ELSE
        CALL CMPTRX(kr,kr,mr,A,Bb,C,B,Tcos,D,W)
      END IF
      DO i = 1, mr
        ii = ip + i
        Q(i,j) = Q(i,jm2) + P(ii) + B(i)
      END DO
    ELSE
      DO i = 1, mr
        B(i) = Q(i,j)
      END DO
      CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
      ip = 0
      ipstor = mr
      IF ( Istag==1 ) THEN
        DO i = 1, mr
          P(i) = B(i)
          Q(i,j) = Q(i,jm2) + 2.*Q(i,jp2) + 3.*B(i)
        END DO
      ELSE
        DO i = 1, mr
          P(i) = B(i)
          B(i) = B(i) + Q(i,N)
        END DO
        Tcos(1) = CMPLX(1.,0.)
        Tcos(2) = CMPLX(0.,0.)
        CALL CMPTRX(1,1,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,j) = Q(i,jm2) + P(i) + B(i)
        END DO
      END IF
    END IF
    lr = kr
    kr = kr + jr
  ELSE
    !
    !     ODD NUMBER OF UNKNOWNS
    !
    IF ( i2r/=1 ) THEN
      DO i = 1, mr
        B(i) = Q(i,j) + .5*(Q(i,jm2)-Q(i,jm1)-Q(i,jm3))
      END DO
      IF ( nrodpr/=0 ) THEN
        DO i = 1, mr
          Q(i,j) = Q(i,j) - Q(i,jm1) + Q(i,jm2)
        END DO
      ELSE
        DO i = 1, mr
          ii = ip + i
          Q(i,j) = Q(i,jm2) + P(ii)
        END DO
        ip = ip - mr
      END IF
      IF ( lr==0 ) THEN
        DO i = 1, mr
          B(i) = fistag*B(i)
        END DO
      ELSE
        CALL CMPCSG(lr,1,0.5,fden,Tcos(kr+1))
      END IF
    ELSE
      DO i = 1, mr
        B(i) = fistag*Q(i,j)
        Q(i,j) = Q(i,jm2)
      END DO
    END IF
    CALL CMPCSG(kr,1,0.5,fden,Tcos)
    CALL CMPTRX(kr,lr,mr,A,Bb,C,B,Tcos,D,W)
    DO i = 1, mr
      Q(i,j) = Q(i,j) + B(i)
    END DO
    kr = kr + i2r
  END IF
  IF ( Mixbnd==2 ) THEN
    nr = nlast/jr
    IF ( nr<=1 ) THEN
      DO i = 1, mr
        B(i) = Q(i,nlast)
      END DO
      GOTO 500
    END IF
  ELSE
    nr = (nlast-1)/jr + 1
    IF ( nr<=3 ) GOTO 200
  END IF
  i2r = jr
  nrodpr = nrod
  GOTO 100
  !
  !      BEGIN SOLUTION
  !
  200  j = 1 + jr
  jm1 = j - i2r
  jp1 = j + i2r
  jm2 = nlast - i2r
  IF ( nr==2 ) THEN
    IF ( N/=2 ) THEN
      !
      !     CASE OF GENERAL N AND NR = 2 .
      !
      DO i = 1, mr
        ii = ip + i
        B3(i) = (0.,0.)
        B(i) = Q(i,1) + 2.*P(ii)
        Q(i,1) = .5*Q(i,1) - Q(i,jm1)
        B2(i) = 2.*(Q(i,1)+Q(i,nlast))
      END DO
      k(1) = kr + jr - 1
      Tcos(k(1)+1) = (-2.,0.)
      k(4) = k(1) + 3 - Istag
      CALL CMPCSG(kr+Istag-2,1,0.0,fnum,Tcos(k(4)))
      k(4) = k(1) + kr + 1
      CALL CMPCSG(jr-1,1,0.0,1.0,Tcos(k(4)))
      CALL C1MERG(Tcos,k(1),kr,k(1)+kr,jr-1,0)
      CALL CMPCSG(kr,1,0.5,fden,Tcos(k(1)+1))
      k(2) = kr
      k(4) = k(1) + k(2) + 1
      CALL CMPCSG(lr,1,0.5,fden,Tcos(k(4)))
      k(3) = lr
      k(4) = 0
      CALL CMPTR3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
      DO i = 1, mr
        B(i) = B(i) + B2(i)
      END DO
      Tcos(1) = (2.,0.)
      CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,1) = Q(i,1) + B(i)
      END DO
    ELSE
      !
      !     CASE  N = 2
      !
      DO i = 1, mr
        B(i) = Q(i,1)
      END DO
      Tcos(1) = (0.,0.)
      CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,1) = B(i)
        B(i) = 2.*(Q(i,2)+B(i))*fistag
      END DO
      Tcos(1) = CMPLX(-fistag,0.)
      Tcos(2) = CMPLX(2.,0.)
      CALL CMPTRX(2,0,mr,A,Bb,C,B,Tcos,D,W)
      DO i = 1, mr
        Q(i,1) = Q(i,1) + B(i)
      END DO
      jr = 1
      i2r = 0
    END IF
    GOTO 400
  ELSE
    IF ( lr==0 ) THEN
      IF ( N/=3 ) THEN
        !
        !     CASE N = 2**P+1
        !
        IF ( Istag/=2 ) THEN
          DO i = 1, mr
            B(i) = Q(i,j) + .5*Q(i,1) - Q(i,jm1) + Q(i,nlast) - Q(i,jm2)
          END DO
          CALL CMPCSG(jr,1,0.5,0.0,Tcos)
          CALL CMPTRX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
          DO i = 1, mr
            Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1)) + B(i)
            B(i) = Q(i,1) + 2.*Q(i,nlast) + 4.*Q(i,j)
          END DO
          jr2 = 2*jr
          CALL CMPCSG(jr,1,0.0,0.0,Tcos)
          DO i = 1, jr
            i1 = jr + i
            i2 = jr + 1 - i
            Tcos(i1) = -Tcos(i2)
          END DO
          CALL CMPTRX(jr2,0,mr,A,Bb,C,B,Tcos,D,W)
          DO i = 1, mr
            Q(i,j) = Q(i,j) + B(i)
            B(i) = Q(i,1) + 2.*Q(i,j)
          END DO
          CALL CMPCSG(jr,1,0.5,0.0,Tcos)
          CALL CMPTRX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
          DO i = 1, mr
            Q(i,1) = .5*Q(i,1) - Q(i,jm1) + B(i)
          END DO
          GOTO 400
        END IF
        !
        !     CASE N = 3.
        !
      ELSEIF ( Istag==2 ) THEN
        !
        !     CASE OF GENERAL N WITH NR = 3 .
        !
        DO i = 1, mr
          B(i) = Q(i,2)
          Q(i,2) = (0.,0.)
          B2(i) = Q(i,3)
          B3(i) = Q(i,1)
        END DO
        jr = 1
        i2r = 0
        j = 2
        GOTO 300
      ELSE
        DO i = 1, mr
          B(i) = Q(i,2)
        END DO
        Tcos(1) = CMPLX(0.,0.)
        CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,2) = B(i)
          B(i) = 4.*B(i) + Q(i,1) + 2.*Q(i,3)
        END DO
        Tcos(1) = CMPLX(-2.,0.)
        Tcos(2) = CMPLX(2.,0.)
        i1 = 2
        i2 = 0
        CALL CMPTRX(i1,i2,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,2) = Q(i,2) + B(i)
          B(i) = Q(i,1) + 2.*Q(i,2)
        END DO
        Tcos(1) = (0.,0.)
        CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,1) = B(i)
        END DO
        jr = 1
        i2r = 0
        GOTO 400
      END IF
    END IF
    DO i = 1, mr
      B(i) = .5*Q(i,1) - Q(i,jm1) + Q(i,j)
    END DO
    IF ( nrod/=0 ) THEN
      DO i = 1, mr
        B(i) = B(i) + Q(i,nlast) - Q(i,jm2)
      END DO
    ELSE
      DO i = 1, mr
        ii = ip + i
        B(i) = B(i) + P(ii)
      END DO
    END IF
    DO i = 1, mr
      t = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
      Q(i,j) = t
      B2(i) = Q(i,nlast) + t
      B3(i) = Q(i,1) + 2.*t
    END DO
  END IF
  300  k(1) = kr + 2*jr - 1
  k(2) = kr + jr
  Tcos(k(1)+1) = (-2.,0.)
  k(4) = k(1) + 3 - Istag
  CALL CMPCSG(k(2)+Istag-2,1,0.0,fnum,Tcos(k(4)))
  k(4) = k(1) + k(2) + 1
  CALL CMPCSG(jr-1,1,0.0,1.0,Tcos(k(4)))
  CALL C1MERG(Tcos,k(1),k(2),k(1)+k(2),jr-1,0)
  k(3) = k(1) + k(2) + lr
  CALL CMPCSG(jr,1,0.5,0.0,Tcos(k(3)+1))
  k(4) = k(3) + jr + 1
  CALL CMPCSG(kr,1,0.5,fden,Tcos(k(4)))
  CALL C1MERG(Tcos,k(3),jr,k(3)+jr,kr,k(1))
  IF ( lr/=0 ) THEN
    CALL CMPCSG(lr,1,0.5,fden,Tcos(k(4)))
    CALL C1MERG(Tcos,k(3),jr,k(3)+jr,lr,k(3)-lr)
    CALL CMPCSG(kr,1,0.5,fden,Tcos(k(4)))
  END IF
  k(3) = kr
  k(4) = kr
  CALL CMPTR3(mr,A,Bb,C,k,B,B2,B3,Tcos,D,W,W2,W3)
  DO i = 1, mr
    B(i) = B(i) + B2(i) + B3(i)
  END DO
  Tcos(1) = (2.,0.)
  CALL CMPTRX(1,0,mr,A,Bb,C,B,Tcos,D,W)
  DO i = 1, mr
    Q(i,j) = Q(i,j) + B(i)
    B(i) = Q(i,1) + 2.*Q(i,j)
  END DO
  CALL CMPCSG(jr,1,0.5,0.0,Tcos)
  CALL CMPTRX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
  IF ( jr/=1 ) THEN
    DO i = 1, mr
      Q(i,1) = .5*Q(i,1) - Q(i,jm1) + B(i)
    END DO
  ELSE
    DO i = 1, mr
      Q(i,1) = B(i)
    END DO
  END IF
  !
  !     START BACK SUBSTITUTION.
  !
  400  j = nlast - jr
  DO i = 1, mr
    B(i) = Q(i,nlast) + Q(i,j)
  END DO
  500  jm2 = nlast - i2r
  IF ( jr==1 ) THEN
    DO i = 1, mr
      Q(i,nlast) = (0.,0.)
    END DO
  ELSEIF ( nrod/=0 ) THEN
    DO i = 1, mr
      Q(i,nlast) = Q(i,nlast) - Q(i,jm2)
    END DO
  ELSE
    DO i = 1, mr
      ii = ip + i
      Q(i,nlast) = P(ii)
    END DO
    ip = ip - mr
  END IF
  CALL CMPCSG(kr,1,0.5,fden,Tcos)
  CALL CMPCSG(lr,1,0.5,fden,Tcos(kr+1))
  IF ( lr==0 ) THEN
    DO i = 1, mr
      B(i) = fistag*B(i)
    END DO
  END IF
  CALL CMPTRX(kr,lr,mr,A,Bb,C,B,Tcos,D,W)
  DO i = 1, mr
    Q(i,nlast) = Q(i,nlast) + B(i)
  END DO
  nlastp = nlast
  DO
    jstep = jr
    jr = i2r
    i2r = i2r/2
    IF ( jr==0 ) THEN
      !
      !     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
      !
      W(1) = CMPLX(REAL(ipstor),0.)
      EXIT
    ELSE
      IF ( Mixbnd==2 ) THEN
        jstart = jr
      ELSE
        jstart = 1 + jr
      END IF
      kr = kr - jr
      IF ( nlast+jr>N ) THEN
        jstop = nlast - jr
      ELSE
        kr = kr - jr
        nlast = nlast + jr
        jstop = nlast - jstep
      END IF
      lr = kr - jr
      CALL CMPCSG(jr,1,0.5,0.0,Tcos)
      DO j = jstart, jstop, jstep
        jm2 = j - jr
        jp2 = j + jr
        IF ( j/=jr ) THEN
          DO i = 1, mr
            B(i) = Q(i,j) + Q(i,jm2) + Q(i,jp2)
          END DO
        ELSE
          DO i = 1, mr
            B(i) = Q(i,j) + Q(i,jp2)
          END DO
        END IF
        IF ( jr/=1 ) THEN
          jm1 = j - i2r
          jp1 = j + i2r
          DO i = 1, mr
            Q(i,j) = .5*(Q(i,j)-Q(i,jm1)-Q(i,jp1))
          END DO
        ELSE
          DO i = 1, mr
            Q(i,j) = (0.,0.)
          END DO
        END IF
        CALL CMPTRX(jr,0,mr,A,Bb,C,B,Tcos,D,W)
        DO i = 1, mr
          Q(i,j) = Q(i,j) + B(i)
        END DO
      END DO
      nrod = 1
      IF ( nlast+i2r<=N ) nrod = 0
      IF ( nlastp/=nlast ) GOTO 400
    END IF
  END DO
END SUBROUTINE CMPOSN
