!** DDNTL
PURE SUBROUTINE DDNTL(Eps,F,FA,Hmax,Hold,Impl,Jtask,Matdim,Maxord,Mint,Miter,&
    Ml,Mu,N,Nde,Save1,T,Uround,USERS,Y,Ywt,H,Mntold,Mtrold,Nfe,Rc,Yh,A,Convrg,&
    El,Fac,Ier,Ipvt,Nq,Nwait,Rh,Rmax,Save2,Tq,Trend,Iswflg,Jstate)
  !> Subroutine DDNTL is called to set parameters on the first call to DDSTP,
  !  on an internal restart, or when the user has altered MINT, MITER, and/or H.
  !***
  ! **Library:**   SLATEC (SDRIVE)
  !***
  ! **Type:**      DOUBLE PRECISION (SDNTL-S, DDNTL-D, CDNTL-C)
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
  !  If the caller has changed MINT, or if JTASK = 0, DDCST is called
  !  to set the coefficients of the method.  If the caller has changed H,
  !  YH must be rescaled.  If H or MINT has been changed, NWAIT is
  !  reset to NQ + 2 to prevent further increases in H for that many
  !  steps.  Also, RC is reset.  RC is the ratio of new to old values of
  !  the coefficient L(0)*H.  If the caller has changed MITER, RC is
  !  set to 0 to force the partials to be updated, if partials are used.
  !
  !***
  ! **Routines called:**  DDCST, DDSCL, DGBFA, DGBSL, DGEFA, DGESL, DNRM2

  !* REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  USE linpack, ONLY : DGBFA, DGEFA
  USE lapack, ONLY : DGBTRS, DGETRS
  !
  INTERFACE
    PURE SUBROUTINE F(N,T,Y,Ydot)
      IMPORT DP
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: T, Y(:)
      REAL(DP), INTENT(OUT) :: Ydot(:)
    END SUBROUTINE F
    PURE SUBROUTINE USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,Iflag)
      IMPORT DP
      INTEGER, INTENT(IN) :: Impl, N, Nde, iflag
      REAL(DP), INTENT(IN) :: T, H, El
      REAL(DP), INTENT(IN) :: Y(N), Yh(N,13), Ywt(N)
      REAL(DP), INTENT(INOUT) :: Save1(N), Save2(N)
    END SUBROUTINE USERS
    PURE SUBROUTINE FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
      IMPORT DP
      INTEGER, INTENT(IN) :: N, Matdim, Ml, Mu, Nde
      REAL(DP), INTENT(IN) :: T, Y(N)
      REAL(DP), INTENT(INOUT) :: A(:,:)
    END SUBROUTINE FA
  END INTERFACE
  INTEGER, INTENT(IN) :: Impl, Iswflg, Jtask, Matdim, Maxord, Mint, Miter, Ml, &
    Mu, N, Nde
  INTEGER, INTENT(INOUT) :: Mntold, Mtrold, Nfe, Nq
  INTEGER, INTENT(OUT) :: Jstate, Nwait
  INTEGER, INTENT(OUT) :: Ipvt(N)
  REAL(DP), INTENT(IN) :: Eps, Hmax, T, Uround
  REAL(DP), INTENT(INOUT) :: H, Hold, Rc, Rh, Rmax, Trend
  REAL(DP), INTENT(INOUT) :: El(13,12), Tq(3,12)
  REAL(DP), INTENT(IN) :: Y(N+1), Ywt(N)
  REAL(DP), INTENT(INOUT) :: A(Matdim,N), Fac(N), Save1(N), Save2(N), Yh(N,13)
  LOGICAL, INTENT(INOUT) :: Convrg
  LOGICAL, INTENT(OUT) :: Ier
  !
  INTEGER :: i, iflag, info
  REAL(DP) :: oldl0, summ
  REAL(DP), PARAMETER :: RMINIT = 10000._DP
  !* FIRST EXECUTABLE STATEMENT  DDNTL
  Ier = .FALSE.
  IF( Jtask>=0 ) THEN
    IF( Jtask==0 ) THEN
      CALL DDCST(Maxord,Mint,Iswflg,El,Tq)
      Rmax = RMINIT
    END IF
    Rc = 0._DP
    Convrg = .FALSE.
    Trend = 1._DP
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
          CALL DGEFA(A,Matdim,N,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL DGETRS('N',N,1,A,Matdim,Ipvt,Save2,N,info)
        ELSEIF( Miter==4 .OR. Miter==5 ) THEN
          CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL DGBFA(A,Matdim,N,Ml,Mu,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL DGBTRS('N',N,Ml,Mu,1,A,Matdim,Ipvt,Save2,N,info)
        END IF
      ELSEIF( Impl==2 ) THEN
        CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
        IF( N==0 ) THEN
          Jstate = 9
          RETURN
        END IF
        DO i = 1, Nde
          IF( A(i,1)==0._DP ) THEN
            Ier = .TRUE.
            RETURN
          ELSE
            Save2(i) = Save2(i)/A(i,1)
          END IF
        END DO
        DO i = Nde + 1, N
          A(i,1) = 0._DP
        END DO
      ELSEIF( Impl==3 ) THEN
        IF( Miter==1 .OR. Miter==2 ) THEN
          CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL DGEFA(A,Matdim,Nde,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL DGETRS('N',Nde,1,A,Matdim,Ipvt,Save2,Nde,info)
        ELSEIF( Miter==4 .OR. Miter==5 ) THEN
          CALL FA(N,T,Y,A(Ml+1:,:),Matdim,Ml,Mu,Nde)
          IF( N==0 ) THEN
            Jstate = 9
            RETURN
          END IF
          CALL DGBFA(A,Matdim,Nde,Ml,Mu,Ipvt,info)
          IF( info/=0 ) THEN
            Ier = .TRUE.
            RETURN
          END IF
          CALL DGBTRS('N',Nde,Ml,Mu,1,A,Matdim,Ipvt,Save2,Nde,info)
        END IF
      END IF
    END IF
    DO i = 1, Nde
      Save1(i) = Save2(i)/MAX(1._DP,Ywt(i))
    END DO
    summ = NORM2(Save1(1:Nde))/SQRT(REAL(Nde, DP))
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
      Rc = 0._DP
      Convrg = .FALSE.
    END IF
    IF( Mint/=Mntold ) THEN
      Mntold = Mint
      oldl0 = El(1,Nq)
      CALL DDCST(Maxord,Mint,Iswflg,El,Tq)
      Rc = Rc*El(1,Nq)/oldl0
      Nwait = Nq + 2
    END IF
    IF( H/=Hold ) THEN
      Nwait = Nq + 2
      Rh = H/Hold
      CALL DDSCL(Hmax,N,Nq,Rmax,Hold,Rc,Rh,Yh)
    END IF
  END IF
  !
END SUBROUTINE DDNTL