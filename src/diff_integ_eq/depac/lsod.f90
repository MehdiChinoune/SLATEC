!** LSOD
SUBROUTINE LSOD(F,Neq,T,Y,Tout,Rtol,Atol,Idid,Ypout,Yh,Yh1,Ewt,Savf,Acor,&
    Wm,Iwm,JAC,Intout,Tstop,Tolfac,Delsgn,Rpar,Ipar)
  !>
  !  Subsidiary to DEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (LSOD-S, DLSOD-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !   DEBDF  merely allocates storage for  LSOD  to relieve the user of
  !   the inconvenience of a long call list.  Consequently  LSOD  is used
  !   as described in the comments for  DEBDF .
  !
  !***
  ! **See also:**  DEBDF
  !***
  ! **Routines called:**  HSTART, INTYD, R1MACH, STOD, VNWRMS, XERMSG
  !***
  ! COMMON BLOCKS    DEBDF1

  !* REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  USE DEBDF1, ONLY : TOLd, H, HMIn, HMXi, X => TN, U => UROund, IQUit, INIt, &
    KSTeps, IBEgin, ITOl, IINteg, ITStop, IJAc, IBAnd, JSTart, KFLag, METh, &
    MITer, MAXord, N, NQ, NST, NFE, NJE
  USE service, ONLY : XERMSG, R1MACH
  INTEGER :: Neq, Idid
  INTEGER :: Ipar(:), Iwm(:)
  REAL :: Delsgn, T, Tolfac, Tout, Tstop
  REAL :: Acor(Neq), Atol(:), Ewt(Neq), Rpar(:), Rtol(:), Savf(Neq), Wm(:), Y(Neq), &
    Yh(Neq,6) , Yh1(6*Neq), Ypout(Neq)
  LOGICAL :: Intout
  EXTERNAL :: F, JAC
  INTEGER :: ltol, natolp, nrtolp, intflg, k, l
  REAL :: absdel, big, del, dt, ha, tol
  CHARACTER(8) :: xern1
  CHARACTER(16) :: xern3, xern4
  !
  !.......................................................................
  !
  !  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
  !  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
  !  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
  !  WORK.
  !
  INTEGER, PARAMETER :: maxnum = 500
  !
  !.......................................................................
  !
  !* FIRST EXECUTABLE STATEMENT  LSOD
  IF ( IBEgin==0 ) THEN
    !
    !        ON THE FIRST CALL, PERFORM INITIALIZATION --
    !        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
    !        FUNCTION ROUTINE R1MACH. THE USER MUST MAKE SURE THAT THE
    !        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
    !
    U = R1MACH(4)
    !                          -- SET ASSOCIATED MACHINE DEPENDENT PARAMETER
    Wm(1) = SQRT(U)
    !                          -- SET TERMINATION FLAG
    IQUit = 0
    !                          -- SET INITIALIZATION INDICATOR
    INIt = 0
    !                          -- SET COUNTER FOR ATTEMPTED STEPS
    KSTeps = 0
    !                          -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
    Intout = .FALSE.
    !                          -- SET START INDICATOR FOR STOD CODE
    JSTart = 0
    !                          -- SET BDF METHOD INDICATOR
    METh = 2
    !                          -- SET MAXIMUM ORDER FOR BDF METHOD
    MAXord = 5
    !                          -- SET ITERATION MATRIX INDICATOR
    !
    IF ( IJAc==0.AND.IBAnd==0 ) MITer = 2
    IF ( IJAc==1.AND.IBAnd==0 ) MITer = 1
    IF ( IJAc==0.AND.IBAnd==1 ) MITer = 5
    IF ( IJAc==1.AND.IBAnd==1 ) MITer = 4
    !
    !                          -- SET OTHER NECESSARY ITEMS IN COMMON BLOCK
    N = Neq
    NST = 0
    NJE = 0
    HMXi = 0.
    NQ = 1
    H = 1.
    !                          -- RESET IBEGIN FOR SUBSEQUENT CALLS
    IBEgin = 1
  END IF
  !
  !.......................................................................
  !
  !      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
  !
  IF ( Neq<1 ) THEN
    WRITE (xern1,'(I8)') Neq
    CALL XERMSG('SLATEC','LSOD',&
      'IN DEBDF, THE NUMBER OF EQUATIONS MUST BE A POSITIVE INTEGER.$$YOU HAVE CALLED THE CODE WITH NEQ = '//xern1,6,1)
    Idid = -33
  END IF
  !
  nrtolp = 0
  natolp = 0
  DO k = 1, Neq
    IF ( nrtolp<=0 ) THEN
      IF ( Rtol(k)<0. ) THEN
        WRITE (xern1,'(I8)') k
        WRITE (xern3,'(1PE15.6)') Rtol(k)
        CALL XERMSG('SLATEC','LSOD','IN DEBDF, THE RELATIVE ERROR TOLERANCES&
         & MUST BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH RTOL('//xern1//') = '&
         //xern3//'$$IN THE CASE OF VECTOR ERROR TOLERANCES,&
         & NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.',7,1)
        Idid = -33
        IF ( natolp>0 ) EXIT
        nrtolp = 1
      ELSEIF ( natolp>0 ) THEN
        GOTO 50
      END IF
    END IF
    !
    IF ( Atol(k)<0. ) THEN
      WRITE (xern1,'(I8)') k
      WRITE (xern3,'(1PE15.6)') Atol(k)
      CALL XERMSG('SLATEC','LSOD','IN DEBDF, THE ABSOLUTE ERROR TOLERANCES MUST&
        & BE NON-NEGATIVE.$$YOU HAVE CALLED THE CODE WITH ATOL('//xern1//') = '&
        //xern3//'$$IN THE CASE OF VECTOR ERROR TOLERANCES,&
        & NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.',8,1)
      Idid = -33
      IF ( nrtolp>0 ) EXIT
      natolp = 1
    END IF
    50  IF ( ITOl==0 ) EXIT
  END DO
  !
  IF ( ITStop==1 ) THEN
    IF ( SIGN(1.,Tout-T)/=SIGN(1.,Tstop-T).OR.ABS(Tout-T)>ABS(Tstop-T) ) THEN
      WRITE (xern3,'(1PE15.6)') Tout
      WRITE (xern4,'(1PE15.6)') Tstop
      CALL XERMSG('SLATEC','LSOD','IN DEBDF, YOU HAVE CALLED THE CODE WITH TOUT = '&
        //xern3//'$$BUT YOU HAVE ALSO TOLD THE CODE NOT TO INTEGRATE PAST THE POINT&
        & TSTOP = '//xern4//' BY SETTING INFO(4) = 1.  THESE INSTRUCTIONS CONFLICT.',14,1)
      Idid = -33
    END IF
  END IF
  !
  !        CHECK SOME CONTINUATION POSSIBILITIES
  !
  IF ( INIt/=0 ) THEN
    IF ( T==Tout ) THEN
      WRITE (xern3,'(1PE15.6)') T
      CALL XERMSG('SLATEC','LSOD',&
        'IN DEBDF, YOU HAVE CALLED THE CODE WITH T = TOUT = '//&
        xern3//'  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',9,1)
      Idid = -33
    END IF
    !
    IF ( T/=TOLd ) THEN
      WRITE (xern3,'(1PE15.6)') TOLd
      WRITE (xern4,'(1PE15.6)') T
      CALL XERMSG('SLATEC','LSOD',&
        'IN DEBDF, YOU HAVE CHANGED THE VALUE OF T FROM '//xern3//' TO '//xern4//&
        '  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',10,1)
      Idid = -33
    END IF
    !
    IF ( INIt/=1 ) THEN
      IF ( Delsgn*(Tout-T)<0. ) THEN
        WRITE (xern3,'(1PE15.6)') Tout
        CALL XERMSG('SLATEC','LSOD',&
          'IN DEBDF, BY CALLING THE CODE WITH TOUT = '//xern3//&
          ' YOU ARE ATTEMPTING TO CHANGE THE DIRECTION OF INTEGRATION.$$THIS IS NOT ALLOWED WITHOUT RESTARTING.',11,1)
        Idid = -33
      END IF
    END IF
  END IF
  !
  IF ( Idid==(-33) ) THEN
    IF ( IQUit/=(-33) ) THEN
      !                       INVALID INPUT DETECTED
      IQUit = -33
      IBEgin = -1
    ELSE
      CALL XERMSG('SLATEC','LSOD','IN DEBDF, INVALID INPUT WAS DETECTED ON&
        & SUCCESSIVE ENTRIES.  IT IS IMPOSSIBLE TO PROCEED BECAUSE YOU HAVE NOT&
        & CORRECTED THE PROBLEM, SO EXECUTION IS BEING TERMINATED.',12,2)
    END IF
    RETURN
  END IF
  !
  !.......................................................................
  !
  !     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
  !     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
  !     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
  !     100*U WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
  !
  DO k = 1, Neq
    IF ( Rtol(k)+Atol(k)<=0. ) THEN
      Rtol(k) = 100.*U
      Idid = -2
    END IF
    IF ( ITOl==0 ) EXIT
  END DO
  !
  IF ( Idid/=(-2) ) THEN
    !
    !     BRANCH ON STATUS OF INITIALIZATION INDICATOR
    !            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
    !                   AND DIRECTION NOT YET SET
    !            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
    !            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
    !
    IF ( INIt==0 ) THEN
      !
      !.......................................................................
      !
      !     MORE INITIALIZATION --
      !                         -- EVALUATE INITIAL DERIVATIVES
      !
      INIt = 1
      CALL F(T,Y,Yh(1,2),Rpar,Ipar)
      NFE = 1
      IF ( T==Tout ) THEN
        Idid = 2
        DO l = 1, Neq
          Ypout(l) = Yh(l,2)
        END DO
        TOLd = T
        RETURN
      END IF
    ELSEIF ( INIt/=1 ) THEN
      GOTO 100
    END IF
    !
    !                         -- COMPUTE INITIAL STEP SIZE
    !                         -- SAVE SIGN OF INTEGRATION DIRECTION
    !                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
    !                                              X AND YH(*) FOR STOD
    !
    ltol = 1
    DO l = 1, Neq
      IF ( ITOl==1 ) ltol = l
      tol = Rtol(ltol)*ABS(Y(l)) + Atol(ltol)
      IF ( tol==0. ) GOTO 200
      Ewt(l) = tol
    END DO
    !
    big = SQRT(R1MACH(2))
    CALL HSTART(F,Neq,T,Tout,Y,Yh(1,2),Ewt,1,U,big,Yh(1,3),Yh(1,4),Yh(1,5),&
      Yh(1,6),Rpar,Ipar,H)
    !
    Delsgn = SIGN(1.0,Tout-T)
    X = T
    DO l = 1, Neq
      Yh(l,1) = Y(l)
      Yh(l,2) = H*Yh(l,2)
    END DO
    INIt = 2
  ELSE
    !                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
    !                                                SMALL POSITIVE VALUE
    IBEgin = -1
    RETURN
  END IF
  !
  !.......................................................................
  !
  !   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
  !   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
  !
  100  del = Tout - T
  absdel = ABS(del)
  !
  !.......................................................................
  !
  !   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
  !
  DO WHILE ( ABS(X-T)<absdel )
    !
    !   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
    !   EXTRAPOLATE AND RETURN
    !
    IF ( ITStop==1 ) THEN
      IF ( ABS(Tstop-X)<100.*U*ABS(X) ) THEN
        dt = Tout - X
        DO l = 1, Neq
          Y(l) = Yh(l,1) + (dt/H)*Yh(l,2)
        END DO
        CALL F(Tout,Y,Ypout,Rpar,Ipar)
        NFE = NFE + 1
        Idid = 3
        T = Tout
        TOLd = T
        RETURN
      END IF
    END IF
    !
    IF ( .NOT.(IINteg==0.OR..NOT.Intout) ) THEN
      !
      !   INTERMEDIATE-OUTPUT MODE
      !
      Idid = 1
      GOTO 300
      !
      !.......................................................................
      !
      !     MONITOR NUMBER OF STEPS ATTEMPTED
      !
    ELSEIF ( KSTeps<=maxnum ) THEN
      !
      !.......................................................................
      !
      !   LIMIT STEP SIZE AND SET WEIGHT VECTOR
      !
      HMIn = 100.*U*ABS(X)
      ha = MAX(ABS(H),HMIn)
      IF ( ITStop==1 ) ha = MIN(ha,ABS(Tstop-X))
      H = SIGN(ha,H)
      ltol = 1
      DO l = 1, Neq
        IF ( ITOl==1 ) ltol = l
        Ewt(l) = Rtol(ltol)*ABS(Yh(l,1)) + Atol(ltol)
        IF ( Ewt(l)<=0.0 ) GOTO 200
      END DO
      Tolfac = U*VNWRMS(Neq,Yh,Ewt)
      IF ( Tolfac<=1. ) THEN
        !
        !.......................................................................
        !
        !     TAKE A STEP
        !
        CALL STOD(Neq,Y,Yh,Neq,Yh1,Ewt,Savf,Acor,Wm,Iwm,F,JAC,Rpar,Ipar)
        !
        JSTart = -2
        Intout = .TRUE.
        IF ( KFLag/=0 ) THEN
          !
          !.......................................................................
          !
          IF ( KFLag==-1 ) THEN
            !
            !                       REPEATED ERROR TEST FAILURES
            Idid = -7
            IBEgin = -1
          ELSE
            !
            !                       REPEATED CORRECTOR CONVERGENCE FAILURES
            Idid = -6
            IBEgin = -1
          END IF
          GOTO 300
        END IF
      ELSE
        !
        !                       TOLERANCES TOO SMALL
        Idid = -2
        Tolfac = 2.*Tolfac
        Rtol(1) = Tolfac*Rtol(1)
        Atol(1) = Tolfac*Atol(1)
        IF ( ITOl/=0 ) THEN
          DO l = 2, Neq
            Rtol(l) = Tolfac*Rtol(l)
            Atol(l) = Tolfac*Atol(l)
          END DO
        END IF
        IBEgin = -1
        GOTO 300
      END IF
    ELSE
      !
      !                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      Idid = -1
      KSTeps = 0
      IBEgin = -1
      GOTO 300
    END IF
  END DO
  CALL INTYD(Tout,0,Yh,Neq,Y,intflg)
  CALL INTYD(Tout,1,Yh,Neq,Ypout,intflg)
  Idid = 3
  IF ( X==Tout ) THEN
    Idid = 2
    Intout = .FALSE.
  END IF
  T = Tout
  TOLd = T
  RETURN
  !
  !                       RELATIVE ERROR CRITERION INAPPROPRIATE
  200  Idid = -3
  IBEgin = -1
  !
  !.......................................................................
  !
  !                       STORE VALUES BEFORE RETURNING TO DEBDF
  300 CONTINUE
  DO l = 1, Neq
    Y(l) = Yh(l,1)
    Ypout(l) = Yh(l,2)/H
  END DO
  T = X
  TOLd = T
  Intout = .FALSE.
END SUBROUTINE LSOD
