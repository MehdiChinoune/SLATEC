!DECK SDAINI
SUBROUTINE SDAINI(X,Y,Yprime,Neq,RES,JAC,H,Wt,Idid,Rpar,Ipar,Phi,Delta,E,&
    Wm,Iwm,Hmin,Uround,Nonneg,Ntemp)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SDAINI
  !***SUBSIDIARY
  !***PURPOSE  Initialization routine for SDASSL.
  !***LIBRARY   SLATEC (DASSL)
  !***TYPE      SINGLE PRECISION (SDAINI-S, DDAINI-D)
  !***AUTHOR  Petzold, Linda R., (LLNL)
  !***DESCRIPTION
  !-----------------------------------------------------------------
  !     SDAINI TAKES ONE STEP OF SIZE H OR SMALLER
  !     WITH THE BACKWARD EULER METHOD, TO
  !     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
  !     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
  !     SOLVE THE CORRECTOR ITERATION.
  !
  !     THE INITIAL GUESS FOR YPRIME IS USED IN THE
  !     PREDICTION, AND IN FORMING THE ITERATION
  !     MATRIX, BUT IS NOT INVOLVED IN THE
  !     ERROR TEST. THIS MAY HAVE TROUBLE
  !     CONVERGING IF THE INITIAL GUESS IS NO
  !     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
  !     NONLINEARLY ON YPRIME.
  !
  !     THE PARAMETERS REPRESENT:
  !     X --         INDEPENDENT VARIABLE
  !     Y --         SOLUTION VECTOR AT X
  !     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
  !     NEQ --       NUMBER OF EQUATIONS
  !     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
  !                  SMALLER THAN H.
  !     WT --        VECTOR OF WEIGHTS FOR ERROR
  !                  CRITERION
  !     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
  !                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
  !                  IDID=-12 -- SDAINI FAILED TO FIND YPRIME
  !     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
  !                  THAT ARE NOT ALTERED BY SDAINI
  !     PHI --       WORK SPACE FOR SDAINI
  !     DELTA,E --   WORK SPACE FOR SDAINI
  !     WM,IWM --    REAL AND INTEGER ARRAYS STORING
  !                  MATRIX INFORMATION
  !
  !-----------------------------------------------------------------
  !***ROUTINES CALLED  SDAJAC, SDANRM, SDASLV
  !***REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !   901030  Minor corrections to declarations.  (FNF)
  !***END PROLOGUE  SDAINI
  !
  INTEGER Neq, Idid, Ipar(*), Iwm(*), Nonneg, Ntemp
  REAL X, Y(*), Yprime(*), H, Wt(*), Rpar(*), Phi(Neq,*), Delta(*), &
    E(*), Wm(*), Hmin, Uround
  EXTERNAL RES, JAC
  !
  EXTERNAL SDAJAC, SDANRM, SDASLV
  REAL SDANRM
  !
  INTEGER i, ier, ires, jcalc, LNJE, LNRE, m, maxit, mjac, ncf, &
    nef, nsf
  REAL cj, damp, delnrm, err, oldnrm, r, rate, s, xold, ynorm
  LOGICAL convgd
  !
  PARAMETER (LNRE=12)
  PARAMETER (LNJE=13)
  !
  DATA maxit/10/, mjac/5/
  DATA damp/0.75E0/
  !
  !
  !---------------------------------------------------
  !     BLOCK 1.
  !     INITIALIZATIONS.
  !---------------------------------------------------
  !
  !***FIRST EXECUTABLE STATEMENT  SDAINI
  Idid = 1
  nef = 0
  ncf = 0
  nsf = 0
  xold = X
  ynorm = SDANRM(Neq,Y,Wt,Rpar,Ipar)
  !
  !     SAVE Y AND YPRIME IN PHI
  DO i = 1, Neq
    Phi(i,1) = Y(i)
    Phi(i,2) = Yprime(i)
  ENDDO
  !
  !
  !----------------------------------------------------
  !     BLOCK 2.
  !     DO ONE BACKWARD EULER STEP.
  !----------------------------------------------------
  !
  !     SET UP FOR START OF CORRECTOR ITERATION
  100  cj = 1.0E0/H
  X = X + H
  !
  !     PREDICT SOLUTION AND DERIVATIVE
  DO i = 1, Neq
    Y(i) = Y(i) + H*Yprime(i)
  ENDDO
  !
  jcalc = -1
  m = 0
  convgd = .TRUE.
  DO
    !
    !
    !     CORRECTOR LOOP.
    Iwm(LNRE) = Iwm(LNRE) + 1
    ires = 0
    !
    CALL RES(X,Y,Yprime,Delta,ires,Rpar,Ipar)
    IF ( ires<0 ) THEN
      !
      !
      !     EXITS FROM CORRECTOR LOOP.
      convgd = .FALSE.
    ELSE
      !
      !
      !     EVALUATE THE ITERATION MATRIX
      IF ( jcalc==-1 ) THEN
        Iwm(LNJE) = Iwm(LNJE) + 1
        jcalc = 0
        CALL SDAJAC(Neq,X,Y,Yprime,Delta,cj,H,ier,Wt,E,Wm,Iwm,RES,ires,&
          Uround,JAC,Rpar,Ipar,Ntemp)
        !
        s = 1000000.E0
        IF ( ires<0 ) THEN
          convgd = .FALSE.
          EXIT
        ELSEIF ( ier/=0 ) THEN
          convgd = .FALSE.
          EXIT
        ELSE
          nsf = 0
        ENDIF
      ENDIF
      !
      !
      !
      !     MULTIPLY RESIDUAL BY DAMPING FACTOR
      DO i = 1, Neq
        Delta(i) = Delta(i)*damp
      ENDDO
      !
      !     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
      !     STORE THE CORRECTION IN DELTA
      !
      CALL SDASLV(Neq,Delta,Wm,Iwm)
      !
      !     UPDATE Y AND YPRIME
      DO i = 1, Neq
        Y(i) = Y(i) - Delta(i)
        Yprime(i) = Yprime(i) - cj*Delta(i)
      ENDDO
      !
      !     TEST FOR CONVERGENCE OF THE ITERATION.
      !
      delnrm = SDANRM(Neq,Delta,Wt,Rpar,Ipar)
      IF ( delnrm>100.E0*Uround*ynorm ) THEN
        !
        IF ( m>0 ) THEN
          !
          rate = (delnrm/oldnrm)**(1.0E0/m)
          IF ( rate>0.90E0 ) THEN
            convgd = .FALSE.
            EXIT
          ELSE
            s = rate/(1.0E0-rate)
          ENDIF
        ELSE
          oldnrm = delnrm
        ENDIF
        !
        IF ( s*delnrm>0.33E0 ) THEN
          !
          !
          !     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
          !     M AND AND TEST WHETHER THE MAXIMUM
          !     NUMBER OF ITERATIONS HAVE BEEN TRIED.
          !     EVERY MJAC ITERATIONS, GET A NEW
          !     ITERATION MATRIX.
          !
          m = m + 1
          IF ( m>=maxit ) THEN
            convgd = .FALSE.
            EXIT
          ELSE
            !
            IF ( (m/mjac)*mjac==m ) jcalc = -1
            CYCLE
          ENDIF
        ENDIF
      ENDIF
      !
      !
      !     THE ITERATION HAS CONVERGED.
      !     CHECK NONNEGATIVITY CONSTRAINTS
      IF ( Nonneg/=0 ) THEN
        DO i = 1, Neq
          Delta(i) = MIN(Y(i),0.0E0)
        ENDDO
        !
        delnrm = SDANRM(Neq,Delta,Wt,Rpar,Ipar)
        IF ( delnrm>0.33E0 ) THEN
          convgd = .FALSE.
        ELSE
          !
          DO i = 1, Neq
            Y(i) = Y(i) - Delta(i)
            Yprime(i) = Yprime(i) - cj*Delta(i)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    EXIT
  ENDDO
  IF ( convgd ) THEN
    !
    !
    !
    !-----------------------------------------------------
    !     BLOCK 3.
    !     THE CORRECTOR ITERATION CONVERGED.
    !     DO ERROR TEST.
    !-----------------------------------------------------
    !
    DO i = 1, Neq
      E(i) = Y(i) - Phi(i,1)
    ENDDO
    err = SDANRM(Neq,E,Wt,Rpar,Ipar)
    !
    IF ( err<=1.0E0 ) RETURN
  ENDIF
  !
  !
  !
  !--------------------------------------------------------
  !     BLOCK 4.
  !     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
  !     AND YPRIME TO THEIR ORIGINAL VALUES.
  !     REDUCE STEPSIZE AND TRY AGAIN, IF
  !     POSSIBLE.
  !---------------------------------------------------------
  !
  X = xold
  DO i = 1, Neq
    Y(i) = Phi(i,1)
    Yprime(i) = Phi(i,2)
  ENDDO
  !
  IF ( convgd ) THEN
    !
    nef = nef + 1
    r = 0.90E0/(2.0E0*err+0.0001E0)
    r = MAX(0.1E0,MIN(0.5E0,r))
    H = H*r
    IF ( ABS(H)>=Hmin.AND.nef<10 ) GOTO 100
  ELSEIF ( ier==0 ) THEN
    IF ( ires>-2 ) THEN
      ncf = ncf + 1
      H = H*0.25E0
      IF ( ncf<10.AND.ABS(H)>=Hmin ) GOTO 100
      Idid = -12
      RETURN
    ELSE
      Idid = -12
      RETURN
    ENDIF
  ELSE
    nsf = nsf + 1
    H = H*0.25E0
    IF ( nsf<3.AND.ABS(H)>=Hmin ) GOTO 100
    Idid = -12
    RETURN
  ENDIF
  Idid = -12
  RETURN
  !
  !-------------END OF SUBROUTINE SDAINI----------------------
END SUBROUTINE SDAINI
