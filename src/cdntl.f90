!*==CDNTL.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK CDNTL
      SUBROUTINE CDNTL(Eps,F,FA,Hmax,Hold,Impl,Jtask,Matdim,Maxord,Mint,Miter,
     &                 Ml,Mu,N,Nde,Save1,T,Uround,USERS,Y,Ywt,H,Mntold,Mtrold,
     &                 Nfe,Rc,Yh,A,Convrg,El,Fac,Ier,Ipvt,Nq,Nwait,Rh,Rmax,
     &                 Save2,Tq,Trend,Iswflg,Jstate)
      IMPLICIT NONE
!*--CDNTL8
!***BEGIN PROLOGUE  CDNTL
!***SUBSIDIARY
!***PURPOSE  Subroutine CDNTL is called to set parameters on the first
!            call to CDSTP, on an internal restart, or when the user has
!            altered MINT, MITER, and/or H.
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDNTL-S, DDNTL-D, CDNTL-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  On the first call, the order is set to 1 and the initial derivatives
!  are calculated.  RMAX is the maximum ratio by which H can be
!  increased in one step.  It is initially RMINIT to compensate
!  for the small initial H, but then is normally equal to RMNORM.
!  If a failure occurs (in corrector convergence or error test), RMAX
!  is set at RMFAIL for the next increase.
!  If the caller has changed MINT, or if JTASK = 0, CDCST is called
!  to set the coefficients of the method.  If the caller has changed H,
!  YH must be rescaled.  If H or MINT has been changed, NWAIT is
!  reset to NQ + 2 to prevent further increases in H for that many
!  steps.  Also, RC is reset.  RC is the ratio of new to old values of
!  the coefficient L(0)*H.  If the caller has changed MITER, RC is
!  set to 0 to force the partials to be updated, if partials are used.
!
!***ROUTINES CALLED  CDCST, CDSCL, CGBFA, CGBSL, CGEFA, CGESL, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDNTL
      INTEGER i , iflag , Impl , info , Iswflg , Jstate , Jtask , Matdim , 
     &        Maxord , Mint , Miter , Ml , Mntold , Mtrold , Mu , N , Nde , 
     &        Nfe , Nq , Nwait
      COMPLEX A(Matdim,*) , Fac(*) , Save1(*) , Save2(*) , Y(*) , Yh(N,*) , 
     &        Ywt(*)
      REAL El(13,12) , Eps , H , Hmax , Hold , oldl0 , Rc , Rh , Rmax , RMINIT , 
     &     SCNRM2 , sum , T , Tq(3,12) , Trend , Uround
      INTEGER Ipvt(*)
      LOGICAL Convrg , Ier
      PARAMETER (RMINIT=10000.E0)
!***FIRST EXECUTABLE STATEMENT  CDNTL
      Ier = .FALSE.
      IF ( Jtask>=0 ) THEN
        IF ( Jtask==0 ) THEN
          CALL CDCST(Maxord,Mint,Iswflg,El,Tq)
          Rmax = RMINIT
        ENDIF
        Rc = 0.E0
        Convrg = .FALSE.
        Trend = 1.E0
        Nq = 1
        Nwait = 3
        CALL F(N,T,Y,Save2)
        IF ( N==0 ) THEN
          Jstate = 6
          RETURN
        ENDIF
        Nfe = Nfe + 1
        IF ( Impl/=0 ) THEN
          IF ( Miter==3 ) THEN
            iflag = 0
            CALL USERS(Y,Yh,Ywt,Save1,Save2,T,H,El,Impl,N,Nde,iflag)
            IF ( iflag==-1 ) THEN
              Ier = .TRUE.
              RETURN
            ENDIF
            IF ( N==0 ) THEN
              Jstate = 10
              RETURN
            ENDIF
          ELSEIF ( Impl==1 ) THEN
            IF ( Miter==1.OR.Miter==2 ) THEN
              CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
              IF ( N==0 ) THEN
                Jstate = 9
                RETURN
              ENDIF
              CALL CGEFA(A,Matdim,N,Ipvt,info)
              IF ( info/=0 ) THEN
                Ier = .TRUE.
                RETURN
              ENDIF
              CALL CGESL(A,Matdim,N,Ipvt,Save2,0)
            ELSEIF ( Miter==4.OR.Miter==5 ) THEN
              CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
              IF ( N==0 ) THEN
                Jstate = 9
                RETURN
              ENDIF
              CALL CGBFA(A,Matdim,N,Ml,Mu,Ipvt,info)
              IF ( info/=0 ) THEN
                Ier = .TRUE.
                RETURN
              ENDIF
              CALL CGBSL(A,Matdim,N,Ml,Mu,Ipvt,Save2,0)
            ENDIF
          ELSEIF ( Impl==2 ) THEN
            CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
            IF ( N==0 ) THEN
              Jstate = 9
              RETURN
            ENDIF
            DO i = 1 , Nde
              IF ( A(i,1)==0.E0 ) THEN
                Ier = .TRUE.
                RETURN
              ELSE
                Save2(i) = Save2(i)/A(i,1)
              ENDIF
            ENDDO
            DO i = Nde + 1 , N
              A(i,1) = 0.E0
            ENDDO
          ELSEIF ( Impl==3 ) THEN
            IF ( Miter==1.OR.Miter==2 ) THEN
              CALL FA(N,T,Y,A,Matdim,Ml,Mu,Nde)
              IF ( N==0 ) THEN
                Jstate = 9
                RETURN
              ENDIF
              CALL CGEFA(A,Matdim,Nde,Ipvt,info)
              IF ( info/=0 ) THEN
                Ier = .TRUE.
                RETURN
              ENDIF
              CALL CGESL(A,Matdim,Nde,Ipvt,Save2,0)
            ELSEIF ( Miter==4.OR.Miter==5 ) THEN
              CALL FA(N,T,Y,A(Ml+1,1),Matdim,Ml,Mu,Nde)
              IF ( N==0 ) THEN
                Jstate = 9
                RETURN
              ENDIF
              CALL CGBFA(A,Matdim,Nde,Ml,Mu,Ipvt,info)
              IF ( info/=0 ) THEN
                Ier = .TRUE.
                RETURN
              ENDIF
              CALL CGBSL(A,Matdim,Nde,Ml,Mu,Ipvt,Save2,0)
            ENDIF
          ENDIF
        ENDIF
        DO i = 1 , Nde
          Save1(i) = Save2(i)/MAX(1.E0,ABS(Ywt(i)))
        ENDDO
        sum = SCNRM2(Nde,Save1,1)/SQRT(REAL(Nde))
        IF ( sum>Eps/ABS(H) ) H = SIGN(Eps/sum,H)
        DO i = 1 , N
          Yh(i,2) = H*Save2(i)
        ENDDO
        IF ( Miter==2.OR.Miter==5.OR.Iswflg==3 ) THEN
          DO i = 1 , N
            Fac(i) = SQRT(Uround)
          ENDDO
        ENDIF
      ELSE
        IF ( Miter/=Mtrold ) THEN
          Mtrold = Miter
          Rc = 0.E0
          Convrg = .FALSE.
        ENDIF
        IF ( Mint/=Mntold ) THEN
          Mntold = Mint
          oldl0 = El(1,Nq)
          CALL CDCST(Maxord,Mint,Iswflg,El,Tq)
          Rc = Rc*El(1,Nq)/oldl0
          Nwait = Nq + 2
        ENDIF
        IF ( H/=Hold ) THEN
          Nwait = Nq + 2
          Rh = H/Hold
          CALL CDSCL(Hmax,N,Nq,Rmax,Hold,Rc,Rh,Yh)
        ENDIF
      ENDIF
      END SUBROUTINE CDNTL
