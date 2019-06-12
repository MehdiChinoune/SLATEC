!** DDPST
SUBROUTINE DDPST(El,F,FA,H,Impl,JACOBN,Matdim,Miter,Ml,Mu,N,Nde,Nq,Save2,&
    T,USERS,Y,Yh,Ywt,Uround,Nfe,Nje,A,Dfdy,Fac,Ier,Ipvt,Save1,Iswflg,Bnd,Jstate)
  !>
  !  Subroutine DDPST evaluates the Jacobian matrix of the right
  !            hand side of the differential equations.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      DOUBLE PRECISION (SDPST-S, DDPST-D, CDPST-C)
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
  ! **Routines called:**  DGBFA, DGEFA, DNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  USE linpack, ONLY : DGBFA, DGEFA
  INTERFACE
    SUBROUTINE F(N,T,Y,Ydot)
      IMPORT DP
      INTEGER :: N
      REAL(DP) :: T, Y(:), Ydot(:)
    END SUBROUTINE F
    SUBROUTINE JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IMPORT DP
      INTEGER :: N, Matdim, Ml, Mu
      REAL(DP) :: T, Y(N), Dfdy(Matdim,N)
    END SUBROUTINE JACOBN
    SUBROUTINE USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      IMPORT DP
      INTEGER :: Impl, N, Nde, iflag
      REAL(DP) :: T, H, El
      REAL(DP) :: Y(N), Yh(N,13), Ywt(N), Save1(N), Save2(N)
    END SUBROUTINE USERS
    SUBROUTINE FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IMPORT DP
      INTEGER :: N, Matdim, Ml, Mu, Nde
      REAL(DP) :: T, Y(N), A(:,:)
    END SUBROUTINE FA
  END INTERFACE
  INTEGER :: Impl, Iswflg, Jstate, Matdim, Miter, Ml, Mu, N, Nde, Nfe, Nje, Nq
  INTEGER :: Ipvt(N)
  REAL(DP) :: Bnd, H, T, Uround
  REAL(DP) :: A(Matdim,N), Dfdy(Matdim,N), El(13,12), Fac(N), Save1(N), Save2(N), &
    Y(N+1), Yh(N,Nq+1), Ywt(N)
  LOGICAL :: Ier
  INTEGER :: i, iflag, imax, info, j, j2, k, mw
  REAL(DP) :: bl, bp, br, dfdymx, diff, dy, facmin, factor, scalee, yj, ys
  REAL(DP), PARAMETER :: FACMAX = 0.5D0, BU = 0.5D0
  !* FIRST EXECUTABLE STATEMENT  DDPST
  Nje = Nje + 1
  Ier = .FALSE.
  IF ( Miter==1.OR.Miter==2 ) THEN
    IF ( Miter==1 ) THEN
      CALL JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      END IF
      IF ( Iswflg==3 ) Bnd = NORM2(Dfdy(1:N,1:N))
      factor = -El(1,Nq)*H
      DO j = 1, N
        DO i = 1, N
          Dfdy(i,j) = factor*Dfdy(i,j)
        END DO
      END DO
    ELSEIF ( Miter==2 ) THEN
      br = Uround**(.875D0)
      bl = Uround**(.75D0)
      bp = Uround**(-.15D0)
      facmin = Uround**(.78D0)
      DO j = 1, N
        ys = MAX(ABS(Ywt(j)),ABS(Y(j)))
        DO
          dy = Fac(j)*ys
          IF ( dy==0.D0 ) THEN
            IF ( Fac(j)<FACMAX ) THEN
              Fac(j) = MIN(100.D0*Fac(j),FACMAX)
              CYCLE
            ELSE
              dy = ys
            END IF
          END IF
          IF ( Nq==1 ) THEN
            dy = SIGN(dy,Save2(j))
          ELSE
            dy = SIGN(dy,Yh(j,3))
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
          factor = -El(1,Nq)*H/dy
          DO i = 1, N
            Dfdy(i,j) = (Save1(i)-Save2(i))*factor
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
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.D0 ) THEN
            scalee = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scalee ) THEN
              Fac(j) = MAX(facmin,Fac(j)*.5D0)
            ELSEIF ( br*scalee<=diff.AND.diff<=bl*scalee ) THEN
              Fac(j) = MIN(Fac(j)*2.D0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scalee ) THEN
              Fac(j) = MIN(bp*Fac(j),FACMAX)
            END IF
          END IF
          EXIT
        END DO
      END DO
      IF ( Iswflg==3 ) Bnd = NORM2(Dfdy(1:N,1:N))/(-El(1,Nq)*H)
      Nfe = Nfe + N
    END IF
    IF ( Impl==0 ) THEN
      DO i = 1, N
        Dfdy(i,i) = Dfdy(i,i) + 1.D0
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
    CALL DGEFA(Dfdy,Matdim,N,Ipvt,info)
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
      br = Uround**(.875D0)
      bl = Uround**(.75D0)
      bp = Uround**(-.15D0)
      facmin = Uround**(.78D0)
      mw = Ml + Mu + 1
      j2 = MIN(mw,N)
      DO j = 1, j2
        DO k = j, N, mw
          ys = MAX(ABS(Ywt(k)),ABS(Y(k)))
          DO
            dy = Fac(k)*ys
            IF ( dy==0.D0 ) THEN
              IF ( Fac(k)<FACMAX ) THEN
                Fac(k) = MIN(100.D0*Fac(k),FACMAX)
                CYCLE
              ELSE
                dy = ys
              END IF
            END IF
            IF ( Nq==1 ) THEN
              dy = SIGN(dy,Save2(k))
            ELSE
              dy = SIGN(dy,Yh(k,3))
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
          Y(k) = Dfdy(mw,k)
          ys = MAX(ABS(Ywt(k)),ABS(Y(k)))
          dy = Fac(k)*ys
          IF ( dy==0.D0 ) dy = ys
          IF ( Nq==1 ) THEN
            dy = SIGN(dy,Save2(k))
          ELSE
            dy = SIGN(dy,Yh(k,3))
          END IF
          dy = (Y(k)+dy) - Y(k)
          factor = -El(1,Nq)*H/dy
          DO i = MAX(Ml+1,mw+1-k), MIN(mw+N-k,mw+Ml)
            Dfdy(i,k) = factor*(Save1(i+k-mw)-Save2(i+k-mw))
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
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.D0 ) THEN
            scalee = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scalee ) THEN
              Fac(j) = MAX(facmin,Fac(j)*.5D0)
            ELSEIF ( br*scalee<=diff.AND.diff<=bl*scalee ) THEN
              Fac(j) = MIN(Fac(j)*2.D0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scalee ) THEN
              Fac(k) = MIN(bp*Fac(k),FACMAX)
            END IF
          END IF
        END DO
      END DO
      Nfe = Nfe + j2
    END IF
    IF ( Iswflg==3 ) THEN
      dfdymx = 0.D0
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          dfdymx = MAX(dfdymx,ABS(Dfdy(i,j)))
        END DO
      END DO
      Bnd = 0.D0
      IF ( dfdymx/=0.D0 ) THEN
        DO j = 1, N
          DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
            Bnd = Bnd + (Dfdy(i,j)/dfdymx)**2
          END DO
        END DO
        Bnd = dfdymx*SQRT(Bnd)/(-El(1,Nq)*H)
      END IF
    END IF
    IF ( Impl==0 ) THEN
      DO j = 1, N
        Dfdy(mw,j) = Dfdy(mw,j) + 1.D0
      END DO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
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
      CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
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
    CALL DGBFA(Dfdy,Matdim,N,Ml,Mu,Ipvt,info)
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
END SUBROUTINE DDPST
