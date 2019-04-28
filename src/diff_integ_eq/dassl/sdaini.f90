!** SDAINI
SUBROUTINE SDAINI(X,Y,Yprime,Neq,RES,JAC,H,Wt,Idid,Rpar,Ipar,Phi,Delta,E,&
    Wm,Iwm,Hmin,Uround,Nonneg,Ntemp)
  !>
  !  Initialization routine for SDASSL.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      SINGLE PRECISION (SDAINI-S, DDAINI-D)
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !***
  ! **Description:**
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
  !***
  ! **Routines called:**  SDAJAC, SDANRM, SDASLV

  !* REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !   901030  Minor corrections to declarations.  (FNF)

  !
  INTEGER Neq, Idid, Ipar(*), Iwm(*), Nonneg, Ntemp
  REAL X, Y(*), Yprime(*), H, Wt(*), Rpar(*), Phi(Neq,*), Delta(*), &
    E(*), Wm(*), Hmin, Uround
  EXTERNAL :: RES, JAC
  !
  INTEGER i, ier, ires, jcalc, m, ncf, nef, nsf
  REAL cj, delnrm, err, oldnrm, r, rate, s, xold, ynorm
  LOGICAL convgd
  !
  INTEGER, PARAMETER :: LNRE = 12
  INTEGER, PARAMETER :: LNJE = 13
  !
  INTEGER, PARAMETER :: maxit = 10, mjac = 5
  REAL, PARAMETER :: damp = 0.75E0
  !
  !
  !---------------------------------------------------
  !     BLOCK 1.
  !     INITIALIZATIONS.
  !---------------------------------------------------
  !
  !* FIRST EXECUTABLE STATEMENT  SDAINI
  Idid = 1
  nef = 0
  ncf = 0
  nsf = 0
  xold = X
  ynorm = SDANRM(Neq,Y,Wt)
  !
  !     SAVE Y AND YPRIME IN PHI
  DO i = 1, Neq
    Phi(i,1) = Y(i)
    Phi(i,2) = Yprime(i)
  END DO
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
  END DO
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
        END IF
      END IF
      !
      !
      !
      !     MULTIPLY RESIDUAL BY DAMPING FACTOR
      DO i = 1, Neq
        Delta(i) = Delta(i)*damp
      END DO
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
      END DO
      !
      !     TEST FOR CONVERGENCE OF THE ITERATION.
      !
      delnrm = SDANRM(Neq,Delta,Wt)
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
          END IF
        ELSE
          oldnrm = delnrm
        END IF
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
          END IF
        END IF
      END IF
      !
      !
      !     THE ITERATION HAS CONVERGED.
      !     CHECK NONNEGATIVITY CONSTRAINTS
      IF ( Nonneg/=0 ) THEN
        DO i = 1, Neq
          Delta(i) = MIN(Y(i),0.0E0)
        END DO
        !
        delnrm = SDANRM(Neq,Delta,Wt)
        IF ( delnrm>0.33E0 ) THEN
          convgd = .FALSE.
        ELSE
          !
          DO i = 1, Neq
            Y(i) = Y(i) - Delta(i)
            Yprime(i) = Yprime(i) - cj*Delta(i)
          END DO
        END IF
      END IF
    END IF
    EXIT
  END DO
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
    END DO
    err = SDANRM(Neq,E,Wt)
    !
    IF ( err<=1.0E0 ) RETURN
  END IF
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
  END DO
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
    END IF
  ELSE
    nsf = nsf + 1
    H = H*0.25E0
    IF ( nsf<3.AND.ABS(H)>=Hmin ) GOTO 100
    Idid = -12
    RETURN
  END IF
  Idid = -12
  RETURN
  !
  !-------------END OF SUBROUTINE SDAINI----------------------
END SUBROUTINE SDAINI
