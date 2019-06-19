!** SDNTL
SUBROUTINE SDNTL(Eps,F,FA,Hmax,Hold,Impl,Jtask,Matdim,Maxord,Mint,Miter,&
    Ml,Mu,N,Nde,Save1,T,Uround,USERS,Y,Ywt,H,Mntold,Mtrold,&
    Nfe,Rc,Yh,A,Convrg,El,Fac,Ier,Ipvt,Nq,Nwait,Rh,Rmax,&
    Save2,Tq,Trend,Iswflg,Jstate)
  !> Subroutine SDNTL is called to set parameters on the first
  !            call to SDSTP, on an internal restart, or when the user has
  !            altered MINT, MITER, and/or H.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      SINGLE PRECISION (SDNTL-S, DDNTL-D, CDNTL-C)
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
  !  On the first call, the order is set to 1 and the initial derivatives
  !  are calculated.  RMAX is the maximum ratio by which H can be
  !  increased in one step.  It is initially RMINIT to compensate
  !  for the small initial H, but then is normally equal to RMNORM.
  !  If a failure occurs (in corrector convergence or error test), RMAX
  !  is set at RMFAIL for the next increase.
  !  If the caller has changed MINT, or if JTASK = 0, SDCST is called
  !  to set the coefficients of the method.  If the caller has changed H,
  !  YH must be rescaled.  If H or MINT has been changed, NWAIT is
  !  reset to NQ + 2 to prevent further increases in H for that many
  !  steps.  Also, RC is reset.  RC is the ratio of new to old values of
  !  the coefficient L(0)*H.  If the caller has changed MITER, RC is
  !  set to 0 to force the partials to be updated, if partials are used.
  !***
  ! **Routines called:**  SDCST, SDSCL, SGBFA, SGBSL, SGEFA, SGESL, SNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  USE linpack, ONLY : SGBFA, SGEFA
  USE lapack, ONLY : SGETRS, SGBTRS
  INTERFACE
    SUBROUTINE F(N,T,Y,Ydot)
      IMPORT SP
      INTEGER :: N
      REAL(SP) :: T, Y(:), Ydot(:)
    END SUBROUTINE F
    SUBROUTINE USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      IMPORT SP
      INTEGER :: Impl, N, Nde, iflag
      REAL(SP) :: T, H, El
      REAL(SP) :: Y(N), Yh(N,13), Ywt(N), Save1(N), Save2(N)
    END SUBROUTINE USERS
    SUBROUTINE FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IMPORT SP
      INTEGER :: N, Matdim, Ml, Mu, Nde
      REAL(SP) :: T, Y(N), A(:,:)
    END SUBROUTINE FA
  END INTERFACE
  INTEGER :: Impl, Iswflg, Jstate, Jtask, Matdim, Maxord, Mint, Miter, Ml, &
    Mntold, Mtrold, Mu, N, Nde, Nfe, Nq, Nwait
  INTEGER :: Ipvt(N)
  REAL(SP) :: Eps, H, Hmax, Hold, Rc, Rh, Rmax, T, Trend, Uround
  REAL(SP) :: A(Matdim,N), El(13,12), Fac(N), Save1(N), Save2(N), Tq(3,12), Y(N+1), &
    Yh(N,13), Ywt(N)
  LOGICAL :: Convrg, Ier
  INTEGER :: i, iflag, info
  REAL(SP) :: oldl0, summ
  REAL(SP), PARAMETER :: RMINIT = 10000.E0
  !* FIRST EXECUTABLE STATEMENT  SDNTL
  Ier = .FALSE.
  IF( Jtask>=0 ) THEN
    IF( Jtask==0 ) THEN
      CALL SDCST(Maxord,Mint,Iswflg,El,Tq)
      Rmax = RMINIT
    END IF
    Rc = 0.E0
    Convrg = .FALSE.
    Trend = 1.E0
    Nq = 1
    Nwait = 3
    CALL F(N,T,Y,Save2)
    IF( N==0 ) THEN
      Jstate = 6
      RETURN
    END IF
    Nfe = Nfe + 1
    IF( Impl/=0 ) THEN
      IF( Miter==3 ) THEN
        iflag = 0
        CALL USERS(Y,Yh,Ywt,Save1,Save2,T,H,El(1,1),Impl,N,Nde,iflag)
        IF( iflag==-1 ) THEN
          Ier = .TRUE.
          RETURN
        END IF
        IF( N==0 ) THEN
          Jstate = 10
          RETURN
        END IF
      ELSEIF( Impl==1 ) THEN
        IF( Miter==1 .OR. Miter==2 ) THEN
          CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL SGEFA(A,Matdim,N,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL SGETRS('N',N,1,A,Matdim,Ipvt,Save2,N,info)
        ELSEIF( Miter==4 .OR. Miter==5 ) THEN
          CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL SGBFA(A,Matdim,N,Ml,Mu,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL SGBTRS('N',N,Ml,Mu,1,A,Matdim,Ipvt,Save2,N,info)
        END IF
      ELSEIF( Impl==2 ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
        DO i = 1, Nde
          IF( A(i,1)==0.E0 ) THEN
            Ier = .TRUE.
            RETURN
          ELSE
            Save2(i) = Save2(i)/A(i,1)
          END IF
        END DO
        DO i = Nde + 1, N
          A(i,1) = 0.E0
        END DO
      ELSEIF( Impl==3 ) THEN
        IF( Miter==1 .OR. Miter==2 ) THEN
          CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL SGEFA(A,Matdim,Nde,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL SGETRS('N',Nde,1,A,Matdim,Ipvt,Save2,Nde,info)
        ELSEIF( Miter==4 .OR. Miter==5 ) THEN
          CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL SGBFA(A,Matdim,Nde,Ml,Mu,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL SGBTRS('N',Nde,Ml,Mu,1,A,Matdim,Ipvt,Save2,Nde,info)
        END IF
      END IF
    END IF
    DO i = 1, Nde
      Save1(i) = Save2(i)/MAX(1.E0,Ywt(i))
    END DO
    summ = NORM2(Save1(1:Nde))/SQRT(REAL(Nde))
    IF( summ>Eps/ABS(H) ) H = SIGN(Eps/summ,H)
    DO i = 1, N
      Yh(i,2) = H*Save2(i)
    END DO
    IF( Miter==2 .OR. Miter==5 .OR. Iswflg==3 ) THEN
      DO i = 1, N
        Fac(i) = SQRT(Uround)
      END DO
    END IF
  ELSE
    IF( Miter/=Mtrold ) THEN
      Mtrold = Miter
      Rc = 0.E0
      Convrg = .FALSE.
    END IF
    IF( Mint/=Mntold ) THEN
      Mntold = Mint
      oldl0 = El(1,Nq)
      CALL SDCST(Maxord,Mint,Iswflg,El,Tq)
      Rc = Rc*El(1,Nq)/oldl0
      Nwait = Nq + 2
    END IF
    IF( H/=Hold ) THEN
      Nwait = Nq + 2
      Rh = H/Hold
      CALL SDSCL(Hmax,N,Nq,Rmax,Hold,Rc,Rh,Yh)
    END IF
  END IF
END SUBROUTINE SDNTL
