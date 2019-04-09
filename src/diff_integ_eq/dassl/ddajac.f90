!** DDAJAC
SUBROUTINE DDAJAC(Neq,X,Y,Yprime,Delta,Cj,H,Ier,Wt,E,Wm,Iwm,RES,Ires,&
    Uround,JAC,Rpar,Ipar,Ntemp)
  IMPLICIT NONE
  !>
  !***
  !  Compute the iteration matrix for DDASSL and form the
  !            LU-decomposition.
  !***
  ! **Library:**   SLATEC (DASSL)
  !***
  ! **Type:**      DOUBLE PRECISION (SDAJAC-S, DDAJAC-D)
  !***
  ! **Author:**  Petzold, Linda R., (LLNL)
  !***
  ! **Description:**
  !-----------------------------------------------------------------------
  !     THIS ROUTINE COMPUTES THE ITERATION MATRIX
  !     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
  !     HERE PD IS COMPUTED BY THE USER-SUPPLIED
  !     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
  !     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
  !     IF IWM(MTYPE)IS 2 OR 5
  !     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
  !     Y        = ARRAY CONTAINING PREDICTED VALUES
  !     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
  !     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
  !                (USED ONLY IF IWM(MTYPE)=2 OR 5)
  !     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
  !     H        = CURRENT STEPSIZE IN INTEGRATION
  !     IER      = VARIABLE WHICH IS .NE. 0
  !                IF ITERATION MATRIX IS SINGULAR,
  !                AND 0 OTHERWISE.
  !     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
  !     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
  !     WM       = REAL WORK SPACE FOR MATRICES. ON
  !                OUTPUT IT CONTAINS THE LU DECOMPOSITION
  !                OF THE ITERATION MATRIX.
  !     IWM      = INTEGER WORK SPACE CONTAINING
  !                MATRIX INFORMATION
  !     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
  !                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
  !     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
  !                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
  !                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
  !                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
  !     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
  !     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
  !                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
  !                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
  !-----------------------------------------------------------------------
  !***
  ! **Routines called:**  DGBFA, DGEFA

  !* REVISION HISTORY  (YYMMDD)
  !   830315  DATE WRITTEN
  !   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
  !   901010  Modified three MAX calls to be all on one line.  (FNF)
  !   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
  !   901026  Added explicit declarations for all variables and minor
  !           cosmetic changes to prologue.  (FNF)
  !   901101  Corrected PURPOSE.  (FNF)

  !
  INTEGER Neq, Ier, Iwm(*), Ires, Ipar(*), Ntemp
  REAL(8) :: X, Y(*), Yprime(*), Delta(*), Cj, H, Wt(*), E(*), Wm(*), Uround, Rpar(*)
  EXTERNAL :: RES, JAC
  !
  EXTERNAL :: DGBFA, DGEFA
  !
  INTEGER i, i1, i2, ii, ipsave, isave, j, k, l, lenpd, mba, mband, meb1, &
    meband, msave, mtype, n, npdm1, nrow
  REAL(8) :: del, delinv, squr, ypsave, ysave
  !
  INTEGER, PARAMETER :: NPD = 1
  INTEGER, PARAMETER :: LML = 1
  INTEGER, PARAMETER :: LMU = 2
  INTEGER, PARAMETER :: LMTYPE = 4
  INTEGER, PARAMETER :: LIPVT = 21
  !
  !* FIRST EXECUTABLE STATEMENT  DDAJAC
  Ier = 0
  npdm1 = NPD - 1
  mtype = Iwm(LMTYPE)
  SELECT CASE (mtype)
    CASE (2)
      !
      !
      !     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
      Ires = 0
      nrow = npdm1
      squr = SQRT(Uround)
      DO i = 1, Neq
        del = squr*MAX(ABS(Y(i)),ABS(H*Yprime(i)),ABS(Wt(i)))
        del = SIGN(del,H*Yprime(i))
        del = (Y(i)+del) - Y(i)
        ysave = Y(i)
        ypsave = Yprime(i)
        Y(i) = Y(i) + del
        Yprime(i) = Yprime(i) + Cj*del
        CALL RES(X,Y,Yprime,E,Ires,Rpar,Ipar)
        IF ( Ires<0 ) RETURN
        delinv = 1.0D0/del
        DO l = 1, Neq
          Wm(nrow+l) = (E(l)-Delta(l))*delinv
        END DO
        nrow = nrow + Neq
        Y(i) = ysave
        Yprime(i) = ypsave
      END DO
    CASE (3)
      !
      !
      !     DUMMY SECTION FOR IWM(MTYPE)=3
      RETURN
    CASE (4)
      !
      !
      !     BANDED USER-SUPPLIED MATRIX
      lenpd = (2*Iwm(LML)+Iwm(LMU)+1)*Neq
      DO i = 1, lenpd
        Wm(npdm1+i) = 0.0D0
      END DO
      CALL JAC(X,Y,Yprime,Wm(NPD),Cj,Rpar,Ipar)
      meband = 2*Iwm(LML) + Iwm(LMU) + 1
      !
      !
      !     DO LU DECOMPOSITION OF BANDED PD
      CALL DGBFA(Wm(NPD),meband,Neq,Iwm(LML),Iwm(LMU),Iwm(LIPVT),Ier)
      RETURN
    CASE (5)
      !
      !
      !     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
      mband = Iwm(LML) + Iwm(LMU) + 1
      mba = MIN(mband,Neq)
      meband = mband + Iwm(LML)
      meb1 = meband - 1
      msave = (Neq/mband) + 1
      isave = Ntemp - 1
      ipsave = isave + msave
      Ires = 0
      squr = SQRT(Uround)
      DO j = 1, mba
        DO n = j, Neq, mband
          k = (n-j)/mband + 1
          Wm(isave+k) = Y(n)
          Wm(ipsave+k) = Yprime(n)
          del = squr*MAX(ABS(Y(n)),ABS(H*Yprime(n)),ABS(Wt(n)))
          del = SIGN(del,H*Yprime(n))
          del = (Y(n)+del) - Y(n)
          Y(n) = Y(n) + del
          Yprime(n) = Yprime(n) + Cj*del
        END DO
        CALL RES(X,Y,Yprime,E,Ires,Rpar,Ipar)
        IF ( Ires<0 ) RETURN
        DO n = j, Neq, mband
          k = (n-j)/mband + 1
          Y(n) = Wm(isave+k)
          Yprime(n) = Wm(ipsave+k)
          del = squr*MAX(ABS(Y(n)),ABS(H*Yprime(n)),ABS(Wt(n)))
          del = SIGN(del,H*Yprime(n))
          del = (Y(n)+del) - Y(n)
          delinv = 1.0D0/del
          i1 = MAX(1,(n-Iwm(LMU)))
          i2 = MIN(Neq,(n+Iwm(LML)))
          ii = n*meb1 - Iwm(LML) + npdm1
          DO i = i1, i2
            Wm(ii+i) = (E(i)-Delta(i))*delinv
          END DO
        END DO
      END DO
      CALL DGBFA(Wm(NPD),meband,Neq,Iwm(LML),Iwm(LMU),Iwm(LIPVT),Ier)
      RETURN
    CASE DEFAULT
      !
      !
      !     DENSE USER-SUPPLIED MATRIX
      lenpd = Neq*Neq
      DO i = 1, lenpd
        Wm(npdm1+i) = 0.0D0
      END DO
      CALL JAC(X,Y,Yprime,Wm(NPD),Cj,Rpar,Ipar)
  END SELECT
  !
  !
  !     DO DENSE-MATRIX LU DECOMPOSITION ON PD
  CALL DGEFA(Wm(NPD),Neq,Neq,Iwm(LIPVT),Ier)
  RETURN
  !------END OF SUBROUTINE DDAJAC------
  RETURN
END SUBROUTINE DDAJAC
