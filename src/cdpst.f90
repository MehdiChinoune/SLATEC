!*==CDPST.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CDPST
SUBROUTINE CDPST(El,F,FA,H,Impl,JACOBN,Matdim,Miter,Ml,Mu,N,Nde,Nq,Save2,&
    T,USERS,Y,Yh,Ywt,Uround,Nfe,Nje,A,Dfdy,Fac,Ier,Ipvt,&
    Save1,Iswflg,Bnd,Jstate)
  IMPLICIT NONE
  !*--CDPST7
  !***BEGIN PROLOGUE  CDPST
  !***SUBSIDIARY
  !***PURPOSE  Subroutine CDPST evaluates the Jacobian matrix of the right
  !            hand side of the differential equations.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      COMPLEX (SDPST-S, DDPST-D, CDPST-C)
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
  !***ROUTINES CALLED  CGBFA, CGEFA, SCNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  CDPST
  INTEGER i, iflag, imax, Impl, info, Iswflg, j, j2, Jstate, k, &
    Matdim, Miter, Ml, Mu, mw, N, Nde, Nfe, Nje, Nq
  COMPLEX A(Matdim,*), cfctr, Dfdy(Matdim,*), dy, Fac(*), Save1(*), &
    Save2(*), Y(*), Yh(N,*), yj, ys, Ywt(*)
  REAL bl, Bnd, bp, br, BU, dfdymx, diff, El(13,12), FACMAX, &
    facmin, factor, H, scale, SCNRM2, T, Uround, zmax, zmin
  INTEGER Ipvt(*)
  LOGICAL Ier
  PARAMETER (FACMAX=.5E0,BU=0.5E0)
  !***FIRST EXECUTABLE STATEMENT  CDPST
  Nje = Nje + 1
  Ier = .FALSE.
  IF ( Miter==1.OR.Miter==2 ) THEN
    IF ( Miter==1 ) THEN
      CALL JACOBN(N,T,Y,Dfdy,Matdim,Ml,Mu)
      IF ( N==0 ) THEN
        Jstate = 8
        RETURN
      ENDIF
      IF ( Iswflg==3 ) Bnd = SCNRM2(N*N,Dfdy,1)
      factor = -El(1,Nq)*H
      DO j = 1, N
        DO i = 1, N
          Dfdy(i,j) = factor*Dfdy(i,j)
        ENDDO
      ENDDO
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
        ENDIF
        DO
          dy = Fac(j)*ys
          IF ( dy==0.E0 ) THEN
            IF ( REAL(Fac(j))<FACMAX ) THEN
              Fac(j) = MIN(100.E0*REAL(Fac(j)),FACMAX)
              CYCLE
            ELSE
              dy = ys
            ENDIF
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
          cfctr = -El(1,Nq)*H/dy
          DO i = 1, N
            Dfdy(i,j) = (Save1(i)-Save2(i))*cfctr
          ENDDO
          !                                                                 Step 1
          diff = ABS(Save2(1)-Save1(1))
          imax = 1
          DO i = 2, N
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            ENDIF
          ENDDO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.E0 ) THEN
            scale = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scale ) THEN
              Fac(j) = MAX(facmin,REAL(Fac(j))*.5E0)
            ELSEIF ( br*scale<=diff.AND.diff<=bl*scale ) THEN
              Fac(j) = MIN(REAL(Fac(j))*2.E0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scale ) THEN
              Fac(j) = MIN(bp*REAL(Fac(j)),FACMAX)
            ENDIF
          ENDIF
          EXIT
        ENDDO
      ENDDO
      IF ( Iswflg==3 ) Bnd = SCNRM2(N*N,Dfdy,1)/(-El(1,Nq)*H)
      Nfe = Nfe + N
    ENDIF
    IF ( Impl==0 ) THEN
      DO i = 1, N
        Dfdy(i,i) = Dfdy(i,i) + 1.E0
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1, N
        DO i = 1, N
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO i = 1, Nde
        Dfdy(i,i) = Dfdy(i,i) + A(i,1)
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1, Nde
        DO i = 1, Nde
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ENDIF
    CALL CGEFA(Dfdy,Matdim,N,Ipvt,info)
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
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = factor*Dfdy(i,j)
        ENDDO
      ENDDO
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
          ENDIF
          DO
            dy = Fac(k)*ys
            IF ( dy==0.E0 ) THEN
              IF ( REAL(Fac(k))<FACMAX ) THEN
                Fac(k) = MIN(100.E0*REAL(Fac(k)),FACMAX)
                CYCLE
              ELSE
                dy = ys
              ENDIF
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
        DO k = j, N, mw
          dy = Y(k) - Dfdy(mw,k)
          Y(k) = Dfdy(mw,k)
          cfctr = -El(1,Nq)*H/dy
          DO i = MAX(Ml+1,mw+1-k), MIN(mw+N-k,mw+Ml)
            Dfdy(i,k) = cfctr*(Save1(i+k-mw)-Save2(i+k-mw))
          ENDDO
          !                                                                 Step 1
          imax = MAX(1,k-Mu)
          diff = ABS(Save2(imax)-Save1(imax))
          DO i = MAX(1,k-Mu) + 1, MIN(k+Ml,N)
            IF ( ABS(Save2(i)-Save1(i))>diff ) THEN
              imax = i
              diff = ABS(Save2(i)-Save1(i))
            ENDIF
          ENDDO
          !                                                                 Step 2
          IF ( MIN(ABS(Save2(imax)),ABS(Save1(imax)))>0.E0 ) THEN
            scale = MAX(ABS(Save2(imax)),ABS(Save1(imax)))
            !                                                                 Step 3
            IF ( diff>BU*scale ) THEN
              Fac(j) = MAX(facmin,REAL(Fac(j))*.5E0)
            ELSEIF ( br*scale<=diff.AND.diff<=bl*scale ) THEN
              Fac(j) = MIN(REAL(Fac(j))*2.E0,FACMAX)
              !                                                                 Step 4
            ELSEIF ( diff<br*scale ) THEN
              Fac(k) = MIN(bp*REAL(Fac(k)),FACMAX)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      Nfe = Nfe + j2
    ENDIF
    IF ( Iswflg==3 ) THEN
      dfdymx = 0.E0
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          zmax = MAX(ABS(REAL(Dfdy(i,j))),ABS(AIMAG(Dfdy(i,j))))
          zmin = MIN(ABS(REAL(Dfdy(i,j))),ABS(AIMAG(Dfdy(i,j))))
          IF ( zmax/=0.E0 )&
            dfdymx = MAX(dfdymx,zmax*SQRT(1.E0+(zmin/zmax)**2))
        ENDDO
      ENDDO
      Bnd = 0.E0
      IF ( dfdymx/=0.E0 ) THEN
        DO j = 1, N
          DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
            Bnd = Bnd + (REAL(Dfdy(i,j))/dfdymx)&
              **2 + (AIMAG(Dfdy(i,j))/dfdymx)**2
          ENDDO
        ENDDO
        Bnd = dfdymx*SQRT(Bnd)/(-El(1,Nq)*H)
      ENDIF
    ENDIF
    IF ( Impl==0 ) THEN
      DO j = 1, N
        Dfdy(mw,j) = Dfdy(mw,j) + 1.E0
      ENDDO
    ELSEIF ( Impl==1 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1, N
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+N-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ELSEIF ( Impl==2 ) THEN
      CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1, Nde
        Dfdy(mw,j) = Dfdy(mw,j) + A(j,1)
      ENDDO
    ELSEIF ( Impl==3 ) THEN
      CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
      IF ( N==0 ) THEN
        Jstate = 9
        RETURN
      ENDIF
      DO j = 1, Nde
        DO i = MAX(Ml+1,mw+1-j), MIN(mw+Nde-j,mw+Ml)
          Dfdy(i,j) = Dfdy(i,j) + A(i,j)
        ENDDO
      ENDDO
    ENDIF
    CALL CGBFA(Dfdy,Matdim,N,Ml,Mu,Ipvt,info)
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
END SUBROUTINE CDPST
