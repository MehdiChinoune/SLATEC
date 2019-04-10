!** CDPST
SUBROUTINE CDPST(El,F,FA,H,Impl,JACOBN,Matdim,Miter,Ml,Mu,N,Nde,Nq,Save2,&
    T,USERS,Y,Yh,Ywt,Uround,Nfe,Nje,A,Dfdy,Fac,Ier,Ipvt,Save1,Iswflg,Bnd,Jstate)
  IMPLICIT NONE
  !>
  !***
  !  Subroutine CDPST evaluates the Jacobian matrix of the right
  !            hand side of the differential equations.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      COMPLEX (SDPST-S, DDPST-D, CDPST-C)
  !***
  ! **Author:**  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***
  ! **Description:**
  !
  !  If MITER is 1, 2, 4, or 5, the matrix
  !  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
  !  decomposition, with the results also stored in DFDY.
  !
  !***
  ! **Routines called:**  CGBFA, CGEFA, SCNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.

  INTEGER i, iflag, imax, Impl, info, Iswflg, j, j2, Jstate, k, &
    Matdim, Miter, Ml, Mu, mw, N, Nde, Nfe, Nje, Nq
  COMPLEX A(Matdim,*), cfctr, Dfdy(Matdim,*), dy, Fac(*), Save1(*), &
    Save2(*), Y(*), Yh(N,*), yj, ys, Ywt(*)
  REAL bl, Bnd, bp, br, dfdymx, diff, El(13,12), facmin, factor, H, scalee, &
    SCNRM2, T, Uround, zmax, zmin
  INTEGER Ipvt(*)
  LOGICAL Ier
  REAL, PARAMETER :: FACMAX = 0.5E0, BU = 0.5E0
  !* FIRST EXECUTABLE STATEMENT  CDPST
  Nje = Nje + 1
  Ier = .FALSE.
  IF ( Miter==1.OR.Miter==2 ) THEN
    IF ( Miter==1 ) THEN
      CALL JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      END IF
      IF ( Iswflg==3 ) Bnd = SCNRM2(N*N,Dfdy,1)
      factor = -El(1,Nq)*H
      DO j = 1, N
        DO i = 1, N
          Dfdy(i,j) = factor*Dfdy(i,j)
        END DO
      END DO
    ELSEIF ( Miter==2 ) THEN
      br = Uround**(.875E0)
      bl = Uround**(.75E0)
      bp = Uround**(-.15E0)
      facmin = Uround**(.78E0)
      DO j = 1, N
        IF ( ABS(Y(j))>ABS(Ywt(j)) ) THEN
          ys = Y(j)
        ELSE
          ys = Ywt(j)
        END IF
        DO
          dy = Fac(j)*ys
          IF ( dy==0.E0 ) THEN
            IF ( REAL(Fac(j))<FACMAX ) THEN
              Fac(j) = MIN(100.E0*REAL(Fac(j)),FACMAX)
              CYCLE
            ELSE
              dy = ys
            END IF
          END IF
          dy = (Y(j)+dy) - Y(j)
          yj = Y(j)
          Y(j) = Y(j) + dy
          CALL F(N,T,Y,Save1)
          IF ( N==0 ) THEN
            Jstate = 6
            RETURN
          END IF
          Y(j) = yj
          cfctr = -El(1,Nq)*H/dy
          DO i = 1, N
            Dfdy(i,j) = (Save1(i)-Save2(i))*cfctr
          END DO
          !                                                                 Step 1
          diff = ABS(Save2(1)-Save1(1))
          imax = 1
          DO i = 2, N
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            END IF
          END DO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.E0 ) THEN
            scalee = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scalee ) THEN
              Fac(j) = MAX(facmin,REAL(Fac(j))*.5E0)
            ELSEIF ( br*scalee<=diff.AND.diff<=bl*scalee ) THEN
              Fac(j) = MIN(REAL(Fac(j))*2.E0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scalee ) THEN
              Fac(j) = MIN(bp*REAL(Fac(j)),FACMAX)
            END IF
          END IF
          EXIT
        END DO
      END DO
      IF ( Iswflg==3 ) Bnd = SCNRM2(N*N,Dfdy,1)/(-El(1,Nq)*H)
      Nfe = Nfe + N
    END IF
    IF ( Impl==0 ) THEN
      DO i = 1, N
        Dfdy(i,i) = Dfdy(i,i) + 1.E0
      END DO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO j = 1, N
        DO i = 1, N
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        END DO
      END DO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO i = 1, Nde
        Dfdy(i,i) = Dfdy(i,i) + A(i,1)
      END DO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO j = 1, Nde
        DO i = 1, Nde
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        END DO
      END DO
    END IF
    CALL CGEFA(Dfdy,Matdim,N,Ipvt,info)
    IF ( info/=0 ) Ier = .TRUE.
  ELSEIF ( Miter==4.OR.Miter==5 ) THEN
    IF ( Miter==4 ) THEN
      CALL JACOBN(N,T,Y,Dfdy(Ml+1,1),Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      END IF
      factor = -El(1,Nq)*H
      mw = Ml + Mu + 1
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = factor*Dfdy(i,j)
        END DO
      END DO
    ELSEIF ( Miter==5 ) THEN
      br = Uround**(.875E0)
      bl = Uround**(.75E0)
      bp = Uround**(-.15E0)
      facmin = Uround**(.78E0)
      mw = Ml + Mu + 1
      j2 = MIN(mw,N)
      DO j = 1, j2
        DO k = j, N, mw
          IF ( ABS(Y(k))>ABS(Ywt(k)) ) THEN
            ys = Y(k)
          ELSE
            ys = Ywt(k)
          END IF
          DO
            dy = Fac(k)*ys
            IF ( dy==0.E0 ) THEN
              IF ( REAL(Fac(k))<FACMAX ) THEN
                Fac(k) = MIN(100.E0*REAL(Fac(k)),FACMAX)
                CYCLE
              ELSE
                dy = ys
              END IF
            END IF
            dy = (Y(k)+dy) - Y(k)
            Dfdy(mw,k) = Y(k)
            Y(k) = Y(k) + dy
            EXIT
          END DO
        END DO
        CALL F(N,T,Y,Save1)
        IF ( N==0 ) THEN
          Jstate = 6
          RETURN
        END IF
        DO k = j, N, mw
          dy = Y(k) - Dfdy(mw,k)
          Y(k) = Dfdy(mw,k)
          cfctr = -El(1,Nq)*H/dy
          DO i = MAX(Ml+1,mw+1-k), MIN(mw+N-k,mw+Ml)
            Dfdy(i,k) = cfctr*(Save1(i+k-mw)-Save2(i+k-mw))
          END DO
          !                                                                 Step 1
          imax = MAX(1,k-Mu)
          diff = ABS(Save2(imax)-Save1(imax))
          DO i = MAX(1,k-Mu) + 1, MIN(k+Ml,N)
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            END IF
          END DO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.E0 ) THEN
            scalee = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scalee ) THEN
              Fac(j) = MAX(facmin,REAL(Fac(j))*.5E0)
            ELSEIF ( br*scalee<=diff.AND.diff<=bl*scalee ) THEN
              Fac(j) = MIN(REAL(Fac(j))*2.E0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scalee ) THEN
              Fac(k) = MIN(bp*REAL(Fac(k)),FACMAX)
            END IF
          END IF
        END DO
      END DO
      Nfe = Nfe + j2
    END IF
    IF ( Iswflg==3 ) THEN
      dfdymx = 0.E0
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          zmax = MAX(ABS(REAL(Dfdy(i,j))),ABS(AIMAG(Dfdy(i,j))))
          zmin = MIN(ABS(REAL(Dfdy(i,j))),ABS(AIMAG(Dfdy(i,j))))
          IF ( zmax/=0.E0 )&
            dfdymx = MAX(dfdymx,zmax*SQRT(1.E0+(zmin/zmax)**2))
        END DO
      END DO
      Bnd = 0.E0
      IF ( dfdymx/=0.E0 ) THEN
        DO j = 1, N
          DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
            Bnd = Bnd + (REAL(Dfdy(i,j))/dfdymx)&
              **2 + (AIMAG(Dfdy(i,j))/dfdymx)**2
          END DO
        END DO
        Bnd = dfdymx*SQRT(Bnd)/(-El(1,Nq)*H)
      END IF
    END IF
    IF ( Impl==0 ) THEN
      DO j = 1, N
        Dfdy(mw,j) = Dfdy(mw,j) + 1.E0
      END DO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        END DO
      END DO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO j = 1, Nde
        Dfdy(mw,j) = Dfdy(mw,j) + A(j,1)
      END DO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      END IF
      DO j = 1, Nde
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+Nde-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        END DO
      END DO
    END IF
    CALL CGBFA(Dfdy,Matdim,N,Ml,Mu,Ipvt,info)
    IF ( info/=0 ) Ier = .TRUE.
  ELSEIF ( Miter==3 ) THEN
    iflag = 1
    CALL USERS(Y,Yh(1,2),Ywt,Save1,Save2,T,H,El(1,Nq),Impl,N,Nde,iflag)
    IF ( iflag==-1 ) THEN
      Ier = .TRUE.
      RETURN
    END IF
    IF ( N==0 ) THEN
      Jstate = 10
      RETURN
    END IF
  END IF
END SUBROUTINE CDPST