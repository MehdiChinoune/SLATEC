!*==DDPST.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DDPST
SUBROUTINE DDPST(El,F,FA,H,Impl,JACOBN,Matdim,Miter,Ml,Mu,N,Nde,Nq,Save2,&
    T,USERS,Y,Yh,Ywt,Uround,Nfe,Nje,A,Dfdy,Fac,Ier,Ipvt,&
    Save1,Iswflg,Bnd,Jstate)
  IMPLICIT NONE
  !*--DDPST7
  !***BEGIN PROLOGUE  DDPST
  !***SUBSIDIARY
  !***PURPOSE  Subroutine DDPST evaluates the Jacobian matrix of the right
  !            hand side of the differential equations.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      DOUBLE PRECISION (SDPST-S, DDPST-D, CDPST-C)
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***DESCRIPTION
  !
  !  If MITER is 1, 2, 4, or 5, the matrix
  !  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
  !  decomposition, with the results also stored in DFDY.
  !
  !***ROUTINES CALLED  DGBFA, DGEFA, DNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  DDPST
  INTEGER i , iflag , imax , Impl , info , Iswflg , j , j2 , Jstate , k , &
    Matdim , Miter , Ml , Mu , mw , N , Nde , Nfe , Nje , Nq
  REAL(8) :: A(Matdim,*) , bl , Bnd , bp , br , BU , Dfdy(Matdim,*) , &
    dfdymx , diff , dy , El(13,12) , Fac(*) , FACMAX , &
    facmin , factor , H , Save1(*) , Save2(*) , scale , &
    DNRM2 , T , Uround , Y(*) , Yh(N,*) , yj , ys , Ywt(*)
  INTEGER Ipvt(*)
  LOGICAL Ier
  PARAMETER (FACMAX=.5D0,BU=0.5D0)
  !***FIRST EXECUTABLE STATEMENT  DDPST
  Nje = Nje + 1
  Ier = .FALSE.
  IF ( Miter==1.OR.Miter==2 ) THEN
    IF ( Miter==1 ) THEN
      CALL JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      ENDIF
      IF ( Iswflg==3 ) Bnd = DNRM2(N*N,Dfdy,1)
      factor = -El(1,Nq)*H
      DO j = 1 , N
        DO i = 1 , N
          Dfdy(i,j) = factor*Dfdy(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Miter==2 ) THEN
      br = Uround**(.875D0)
      bl = Uround**(.75D0)
      bp = Uround**(-.15D0)
      facmin = Uround**(.78D0)
      DO j = 1 , N
        ys = MAX(ABS(Ywt(j)),ABS(Y(j)))
        DO
          dy = Fac(j)*ys
          IF ( dy==0.D0 ) THEN
            IF ( Fac(j)<FACMAX ) THEN
              Fac(j) = MIN(100.D0*Fac(j),FACMAX)
              CYCLE
            ELSE
              dy = ys
            ENDIF
          ENDIF
          IF ( Nq==1 ) THEN
            dy = SIGN(dy,Save2(j))
          ELSE
            dy = SIGN(dy,Yh(j,3))
          ENDIF
          dy = (Y(j)+dy) - Y(j)
          yj = Y(j)
          Y(j) = Y(j) + dy
          CALL F(N,T,Y,Save1)
          IF ( N==0 ) THEN
            Jstate = 6
            RETURN
          ENDIF
          Y(j) = yj
          factor = -El(1,Nq)*H/dy
          DO i = 1 , N
            Dfdy(i,j) = (Save1(i)-Save2(i))*factor
          ENDDO
          !                                                                 Step 1
          diff = ABS(Save2(1)-Save1(1))
          imax = 1
          DO i = 2 , N
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            ENDIF
          ENDDO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.D0 ) THEN
            scale = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scale ) THEN
              Fac(j) = MAX(facmin,Fac(j)*.5D0)
            ELSEIF ( br*scale<=diff.AND.diff<=bl*scale ) THEN
              Fac(j) = MIN(Fac(j)*2.D0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scale ) THEN
              Fac(j) = MIN(bp*Fac(j),FACMAX)
            ENDIF
          ENDIF
          EXIT
        ENDDO
      ENDDO
      IF ( Iswflg==3 ) Bnd = DNRM2(N*N,Dfdy,1)/(-El(1,Nq)*H)
      Nfe = Nfe + N
    ENDIF
    IF ( Impl==0 ) THEN
      DO i = 1 , N
        Dfdy(i,i) = Dfdy(i,i) + 1.D0
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1 , N
        DO i = 1 , N
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO i = 1 , Nde
        Dfdy(i,i) = Dfdy(i,i) + A(i,1)
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1 , Nde
        DO i = 1 , Nde
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ENDIF
    CALL DGEFA(Dfdy,Matdim,N,Ipvt,info)
    IF ( info/=0 ) Ier = .TRUE.
  ELSEIF ( Miter==4.OR.Miter==5 ) THEN
    IF ( Miter==4 ) THEN
      CALL JACOBN(N,T,Y,Dfdy(Ml+1,1),Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      ENDIF
      factor = -El(1,Nq)*H
      mw = Ml + Mu + 1
      DO j = 1 , N
        DO i = MAX(Ml+1,mw+1-j) , MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = factor*Dfdy(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Miter==5 ) THEN
      br = Uround**(.875D0)
      bl = Uround**(.75D0)
      bp = Uround**(-.15D0)
      facmin = Uround**(.78D0)
      mw = Ml + Mu + 1
      j2 = MIN(mw,N)
      DO j = 1 , j2
        DO k = j , N , mw
          ys = MAX(ABS(Ywt(k)),ABS(Y(k)))
          DO
            dy = Fac(k)*ys
            IF ( dy==0.D0 ) THEN
              IF ( Fac(k)<FACMAX ) THEN
                Fac(k) = MIN(100.D0*Fac(k),FACMAX)
                CYCLE
              ELSE
                dy = ys
              ENDIF
            ENDIF
            IF ( Nq==1 ) THEN
              dy = SIGN(dy,Save2(k))
            ELSE
              dy = SIGN(dy,Yh(k,3))
            ENDIF
            dy = (Y(k)+dy) - Y(k)
            Dfdy(mw,k) = Y(k)
            Y(k) = Y(k) + dy
            EXIT
          ENDDO
        ENDDO
        CALL F(N,T,Y,Save1)
        IF ( N==0 ) THEN
          Jstate = 6
          RETURN
        ENDIF
        DO k = j , N , mw
          Y(k) = Dfdy(mw,k)
          ys = MAX(ABS(Ywt(k)),ABS(Y(k)))
          dy = Fac(k)*ys
          IF ( dy==0.D0 ) dy = ys
          IF ( Nq==1 ) THEN
            dy = SIGN(dy,Save2(k))
          ELSE
            dy = SIGN(dy,Yh(k,3))
          ENDIF
          dy = (Y(k)+dy) - Y(k)
          factor = -El(1,Nq)*H/dy
          DO i = MAX(Ml+1,mw+1-k) , MIN(mw+N-k,mw+Ml)
            Dfdy(i,k) = factor*(Save1(i+k-mw)-Save2(i+k-mw))
          ENDDO
          !                                                                 Step 1
          imax = MAX(1,k-Mu)
          diff = ABS(Save2(imax)-Save1(imax))
          DO i = MAX(1,k-Mu) + 1 , MIN(k+Ml,N)
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            ENDIF
          ENDDO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.D0 ) THEN
            scale = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scale ) THEN
              Fac(j) = MAX(facmin,Fac(j)*.5D0)
            ELSEIF ( br*scale<=diff.AND.diff<=bl*scale ) THEN
              Fac(j) = MIN(Fac(j)*2.D0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scale ) THEN
              Fac(k) = MIN(bp*Fac(k),FACMAX)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      Nfe = Nfe + j2
    ENDIF
    IF ( Iswflg==3 ) THEN
      dfdymx = 0.D0
      DO j = 1 , N
        DO i = MAX(Ml+1,mw+1-j) , MIN(mw+N-j,mw+Ml)
          dfdymx = MAX(dfdymx,ABS(Dfdy(i,j)))
        ENDDO
      ENDDO
      Bnd = 0.D0
      IF ( dfdymx/=0.D0 ) THEN
        DO j = 1 , N
          DO i = MAX(Ml+1,mw+1-j) , MIN(mw+N-j,mw+Ml)
            Bnd = Bnd + (Dfdy(i,j)/dfdymx)**2
          ENDDO
        ENDDO
        Bnd = dfdymx*SQRT(Bnd)/(-El(1,Nq)*H)
      ENDIF
    ENDIF
    IF ( Impl==0 ) THEN
      DO j = 1 , N
        Dfdy(mw,j) = Dfdy(mw,j) + 1.D0
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1 , N
        DO i = MAX(Ml+1,mw+1-j) , MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1 , Nde
        Dfdy(mw,j) = Dfdy(mw,j) + A(j,1)
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1 , Nde
        DO i = MAX(Ml+1,mw+1-j) , MIN(mw+Nde-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ENDIF
    CALL DGBFA(Dfdy,Matdim,N,Ml,Mu,Ipvt,info)
    IF ( info/=0 ) Ier = .TRUE.
  ELSEIF ( Miter==3 ) THEN
    iflag = 1
    CALL USERS(Y,Yh(1,2),Ywt,Save1,Save2,T,H,El(1,Nq),Impl,N,Nde,iflag)
    IF ( iflag==-1 ) THEN
      Ier = .TRUE.
      RETURN
    ENDIF
    IF ( N==0 ) THEN
      Jstate = 10
      RETURN
    ENDIF
  ENDIF
END SUBROUTINE DDPST
