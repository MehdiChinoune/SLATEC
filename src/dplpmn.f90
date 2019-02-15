!DECK DPLPMN
SUBROUTINE DPLPMN(DUSRMT,Mrelas,Nvars,Costs,Prgopt,Dattrv,Bl,Bu,Ind,Info,&
    Primal,Duals,Amat,Csc,Colnrm,Erd,Erp,Basmat,Wr,Rz,Rg,&
    Rprim,Rhs,Ww,Lmx,Lbm,Ibasis,Ibb,Imat,Ibrc,Ipr,Iwr)
  IMPLICIT NONE
  REAL DUSRMT
  INTEGER i, ibas, idg, ienter, ileave, Info, iopt, ipage, ipagef ,&
    iplace, isave, itbrc, itlp, j, jstrt, k, key, kprint ,&
    Lbm, LCOl
  INTEGER LENl, LENu, Lmx, LP, lpg, lpr, lpr1, lprg, LROw, Mrelas ,&
    mxitlp, n20046, n20058, n20080, n20098, n20119, n20172 ,&
    n20206, n20247, n20252
  INTEGER n20271, n20276, n20283, n20290, NCP, nerr, np, nparm ,&
    npp, npr004, npr005, npr006, npr007, npr008, npr009 ,&
    npr010, npr011, npr012, npr013, npr014
  INTEGER npr015, nredc, ntries, Nvars, nx0066, nx0091, nx0106
  !***BEGIN PROLOGUE  DPLPMN
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SPLPMN-S, DPLPMN-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     MARVEL OPTION(S).. OUTPUT=YES/NO TO ELIMINATE PRINTED OUTPUT.
  !     THIS DOES NOT APPLY TO THE CALLS TO THE ERROR PROCESSOR.
  !
  !     MAIN SUBROUTINE FOR DSPLP PACKAGE.
  !
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  DASUM, DCOPY, DDOT, DPINCW, DPINIT, DPINTM, DPLPCE,
  !                    DPLPDM, DPLPFE, DPLPFL, DPLPMU, DPLPUP, DPNNZR,
  !                    DPOPT, DPRWPG, DVOUT, IVOUT, LA05BD, SCLOSM, XERMSG
  !***COMMON BLOCKS    LA05DD
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !***END PROLOGUE  DPLPMN
  REAL(8) :: abig, aij, Amat(*), anorm, asmall, Basmat(*) ,&
    Bl(*), Bu(*), Colnrm(*), Costs(*), costsc, Csc(*) ,&
    Dattrv(*), dirnrm, Duals(*), dulnrm, eps, tune ,&
    Erd(*), erdnrm, Erp(*), factor, gg, one, Prgopt(*)&
    , Primal(*), resnrm, Rg(*), Rhs(*), rhsnrm, ropt(07)&
    , Rprim(*), rprnrm, Rz(*), rzj, scalr, scosts ,&
    size, SMAll, theta, tolls, upbnd, uu, Wr(*) ,&
    Ww(*), xlamda, xval, zero, rdum(01), tolabs
  REAL(8) :: DDOT, DASUM
  !
  INTEGER Ibasis(*), Ibb(*), Ibrc(Lbm,2), Imat(*), Ind(*), Ipr(*) ,&
    Iwr(*), intopt(08), idum(01)
  !
  !     ARRAY LOCAL VARIABLES
  !     NAME(LENGTH)          DESCRIPTION
  !
  !     COSTS(NVARS)          COST COEFFICIENTS
  !     PRGOPT( )             OPTION VECTOR
  !     DATTRV( )             DATA TRANSFER VECTOR
  !     PRIMAL(NVARS+MRELAS)  AS OUTPUT IT IS PRIMAL SOLUTION OF LP.
  !                           INTERNALLY, THE FIRST NVARS POSITIONS HOLD
  !                           THE COLUMN CHECK SUMS.  THE NEXT MRELAS
  !                           POSITIONS HOLD THE CLASSIFICATION FOR THE
  !                           BASIC VARIABLES  -1 VIOLATES LOWER
  !                           BOUND, 0 FEASIBLE, +1 VIOLATES UPPER BOUND
  !     DUALS(MRELAS+NVARS)   DUAL SOLUTION. INTERNALLY HOLDS R.H. SIDE
  !                           AS FIRST MRELAS ENTRIES.
  !     AMAT(LMX)             SPARSE FORM OF DATA MATRIX
  !     IMAT(LMX)             SPARSE FORM OF DATA MATRIX
  !     BL(NVARS+MRELAS)      LOWER BOUNDS FOR VARIABLES
  !     BU(NVARS+MRELAS)      UPPER BOUNDS FOR VARIABLES
  !     IND(NVARS+MRELAS)     INDICATOR FOR VARIABLES
  !     CSC(NVARS)            COLUMN SCALING
  !     IBASIS(NVARS+MRELAS)  COLS. 1-MRELAS ARE BASIC, REST ARE NON-BASIC
  !     IBB(NVARS+MRELAS)     INDICATOR FOR NON-BASIC VARS., POLARITY OF
  !                           VARS., AND POTENTIALLY INFINITE VARS.
  !                           IF IBB(J).LT.0, VARIABLE J IS BASIC
  !                           IF IBB(J).GT.0, VARIABLE J IS NON-BASIC
  !                           IF IBB(J).EQ.0, VARIABLE J HAS TO BE IGNORED
  !                           BECAUSE IT WOULD CAUSE UNBOUNDED SOLN.
  !                           WHEN MOD(IBB(J),2).EQ.0, VARIABLE IS AT ITS
  !                           UPPER BOUND, OTHERWISE IT IS AT ITS LOWER
  !                           BOUND
  !     COLNRM(NVARS)         NORM OF COLUMNS
  !     ERD(MRELAS)           ERRORS IN DUAL VARIABLES
  !     ERP(MRELAS)           ERRORS IN PRIMAL VARIABLES
  !     BASMAT(LBM)           BASIS MATRIX FOR HARWELL SPARSE CODE
  !     IBRC(LBM,2)           ROW AND COLUMN POINTERS FOR BASMAT(*)
  !     IPR(2*MRELAS)         WORK ARRAY FOR HARWELL SPARSE CODE
  !     IWR(8*MRELAS)         WORK ARRAY FOR HARWELL SPARSE CODE
  !     WR(MRELAS)            WORK ARRAY FOR HARWELL SPARSE CODE
  !     RZ(NVARS+MRELAS)      REDUCED COSTS
  !     RPRIM(MRELAS)         INTERNAL PRIMAL SOLUTION
  !     RG(NVARS+MRELAS)      COLUMN WEIGHTS
  !     WW(MRELAS)            WORK ARRAY
  !     RHS(MRELAS)           HOLDS TRANSLATED RIGHT HAND SIDE
  !
  !     SCALAR LOCAL VARIABLES
  !     NAME       TYPE         DESCRIPTION
  !
  !     LMX        INTEGER      LENGTH OF AMAT(*)
  !     LPG        INTEGER      LENGTH OF PAGE FOR AMAT(*)
  !     EPS        DOUBLE       MACHINE PRECISION
  !     TUNE       DOUBLE       PARAMETER TO SCALE ERROR ESTIMATES
  !     TOLLS      DOUBLE       RELATIVE TOLERANCE FOR SMALL RESIDUALS
  !     TOLABS     DOUBLE       ABSOLUTE TOLERANCE FOR SMALL RESIDUALS.
  !                             USED IF RELATIVE ERROR TEST FAILS.
  !                             IN CONSTRAINT EQUATIONS
  !     FACTOR     DOUBLE      .01--DETERMINES IF BASIS IS SINGULAR
  !                             OR COMPONENT IS FEASIBLE.  MAY NEED TO
  !                             BE INCREASED TO 1.D0 ON SHORT WORD
  !                             LENGTH MACHINES.
  !     ASMALL     DOUBLE       LOWER BOUND FOR NON-ZERO MAGN. IN AMAT(*)
  !     ABIG       DOUBLE       UPPER BOUND FOR NON-ZERO MAGN. IN AMAT(*)
  !     MXITLP     INTEGER      MAXIMUM NUMBER OF ITERATIONS FOR LP
  !     ITLP       INTEGER      ITERATION COUNTER FOR TOTAL LP ITERS
  !     COSTSC     DOUBLE       COSTS(*) SCALING
  !     SCOSTS     DOUBLE       TEMP LOC. FOR COSTSC.
  !     XLAMDA     DOUBLE       WEIGHT PARAMETER FOR PEN. METHOD.
  !     ANORM      DOUBLE       NORM OF DATA MATRIX AMAT(*)
  !     RPRNRM     DOUBLE       NORM OF THE SOLUTION
  !     DULNRM     DOUBLE       NORM OF THE DUALS
  !     ERDNRM     DOUBLE       NORM OF ERROR IN DUAL VARIABLES
  !     DIRNRM     DOUBLE       NORM OF THE DIRECTION VECTOR
  !     RHSNRM     DOUBLE       NORM OF TRANSLATED RIGHT HAND SIDE VECTOR
  !     RESNRM     DOUBLE       NORM OF RESIDUAL VECTOR FOR CHECKING
  !                             FEASIBILITY
  !     NZBM       INTEGER      NUMBER OF NON-ZEROS IN BASMAT(*)
  !     LBM        INTEGER      LENGTH OF BASMAT(*)
  !     SMALL      DOUBLE       EPS*ANORM  USED IN HARWELL SPARSE CODE
  !     LP         INTEGER      USED IN HARWELL LA05*() PACK AS OUTPUT
  !                             FILE NUMBER. SET=I1MACH(4) NOW.
  !     UU         DOUBLE       0.1--USED IN HARWELL SPARSE CODE
  !                             FOR RELATIVE PIVOTING TOLERANCE.
  !     GG         DOUBLE       OUTPUT INFO FLAG IN HARWELL SPARSE CODE
  !     IPLACE     INTEGER      INTEGER USED BY SPARSE MATRIX CODES
  !     IENTER     INTEGER      NEXT COLUMN TO ENTER BASIS
  !     NREDC      INTEGER      NO. OF FULL REDECOMPOSITIONS
  !     KPRINT     INTEGER      LEVEL OF OUTPUT, =0-3
  !     IDG        INTEGER      FORMAT AND PRECISION OF OUTPUT
  !     ITBRC      INTEGER      NO. OF ITERS. BETWEEN RECALCULATING
  !                             THE ERROR IN THE PRIMAL SOLUTION.
  !     NPP        INTEGER      NO. OF NEGATIVE REDUCED COSTS REQUIRED
  !                             IN PARTIAL PRICING
  !     JSTRT      INTEGER      STARTING PLACE FOR PARTIAL PRICING.
  !
  LOGICAL colscp, savedt, contin, cstscp, unbnd, feas, finite ,&
    found, minprb, redbas, singlr, sizeup, stpedg, trans ,&
    usrbas, zerolv, lopt(08)
  CHARACTER(8) :: xern1, xern2
  EQUIVALENCE (contin,lopt(1))
  EQUIVALENCE (usrbas,lopt(2))
  EQUIVALENCE (sizeup,lopt(3))
  EQUIVALENCE (savedt,lopt(4))
  EQUIVALENCE (colscp,lopt(5))
  EQUIVALENCE (cstscp,lopt(6))
  EQUIVALENCE (minprb,lopt(7))
  EQUIVALENCE (stpedg,lopt(8))
  EQUIVALENCE (idg,intopt(1))
  EQUIVALENCE (ipagef,intopt(2))
  EQUIVALENCE (isave,intopt(3))
  EQUIVALENCE (mxitlp,intopt(4))
  EQUIVALENCE (kprint,intopt(5))
  EQUIVALENCE (itbrc,intopt(6))
  EQUIVALENCE (npp,intopt(7))
  EQUIVALENCE (lprg,intopt(8))
  EQUIVALENCE (eps,ropt(1))
  EQUIVALENCE (asmall,ropt(2))
  EQUIVALENCE (abig,ropt(3))
  EQUIVALENCE (costsc,ropt(4))
  EQUIVALENCE (tolls,ropt(5))
  EQUIVALENCE (tune,ropt(6))
  EQUIVALENCE (tolabs,ropt(7))
  !
  !     COMMON BLOCK USED BY LA05 () PACKAGE..
  COMMON /LA05DD/ SMAll, LP, LENl, LENu, NCP, LROw, LCOl
  EXTERNAL DUSRMT
  !
  !     SET LP=0 SO NO ERROR MESSAGES WILL PRINT WITHIN LA05 () PACKAGE.
  !***FIRST EXECUTABLE STATEMENT  DPLPMN
  LP = 0
  !
  !     THE VALUES ZERO AND ONE.
  zero = 0.D0
  one = 1.D0
  factor = 0.01D0
  lpg = Lmx - (Nvars+4)
  iopt = 1
  Info = 0
  unbnd = .FALSE.
  jstrt = 1
  !
  !     PROCESS USER OPTIONS IN PRGOPT(*).
  !     CHECK THAT ANY USER-GIVEN CHANGES ARE WELL-DEFINED.
  CALL DPOPT(Prgopt,Mrelas,Nvars,Info,Csc,Ibasis,ropt,intopt,lopt)
  IF ( Info<0 ) GOTO 4600
  IF ( .NOT.(contin) ) THEN
    !
    !     INITIALIZE SPARSE DATA MATRIX, AMAT(*) AND IMAT(*).
    CALL DPINTM(Mrelas,Nvars,Amat,Imat,Lmx,ipagef)
  ELSE
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PROCEDURE (RETRIEVE SAVED DATA FROM FILE ISAVE)
    lpr = Nvars + 4
    REWIND isave
    READ (isave) (Amat(i),i=1,lpr), (Imat(i),i=1,lpr)
    key = 2
    ipage = 1
    GOTO 2500
  ENDIF
  !
  !     UPDATE MATRIX DATA AND CHECK BOUNDS FOR CONSISTENCY.
  100  CALL DPLPUP(DUSRMT,Mrelas,Nvars,Prgopt,Dattrv,Bl,Bu,Ind,Info,Amat,Imat,&
    sizeup,asmall,abig)
  IF ( Info<0 ) GOTO 4600
  !
  !++  CODE FOR OUTPUT=YES IS ACTIVE
  IF ( kprint>=1 ) THEN
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !++  CODE FOR OUTPUT=YES IS ACTIVE
    !     PROCEDURE (PRINT PROLOGUE)
    idum(1) = Mrelas
    CALL IVOUT(1,idum,'(''1NUM. OF DEPENDENT VARS., MRELAS'')',idg)
    idum(1) = Nvars
    CALL IVOUT(1,idum,'('' NUM. OF INDEPENDENT VARS., NVARS'')',idg)
    CALL IVOUT(1,idum,'('' DIMENSION OF COSTS(*)='')',idg)
    idum(1) = Nvars + Mrelas
    CALL IVOUT(1,idum,&
      '('' DIMENSIONS OF BL(*),BU(*),IND(*)''/'' PRIMAL(*),DUALS(*) ='')',idg)
    CALL IVOUT(1,idum,'('' DIMENSION OF IBASIS(*)='')',idg)
    idum(1) = lprg + 1
    CALL IVOUT(1,idum,'('' DIMENSION OF PRGOPT(*)='')',idg)
    CALL IVOUT(0,idum,'('' 1-NVARS=INDEPENDENT VARIABLE INDICES.''/&
      &'' (NVARS+1)-(NVARS+MRELAS)=DEPENDENT VARIABLE INDICES.''/&
      &'' CONSTRAINT INDICATORS ARE 1-4 AND MEAN'')',idg)
    CALL IVOUT(0,idum,'('' 1=VARIABLE HAS ONLY LOWER BOUND.''/&
      &'' 2=VARIABLE HAS ONLY UPPER BOUND.''/&
      &'' 3=VARIABLE HAS BOTH BOUNDS.''/&
      &'' 4=VARIABLE HAS NO BOUNDS, IT IS FREE.'')',idg)
    CALL DVOUT(Nvars,Costs,'('' ARRAY OF COSTS'')',idg)
    CALL IVOUT(Nvars+Mrelas,Ind,'('' CONSTRAINT INDICATORS'')',idg)
    CALL DVOUT(Nvars+Mrelas,Bl,&
      '('' LOWER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENTRIES.)'')'&
      ,idg)
    CALL DVOUT(Nvars+Mrelas,Bu,&
      '('' UPPER BOUNDS FOR VARIABLES  (IGNORE UNUSED ENTRIES.)'')',idg)
    IF ( kprint>=2 ) THEN
      CALL IVOUT(0,idum,&
        '(''0NON-BASIC INDICES THAT ARE NEGATIVE SHOW VARIABLES EXCHANGED AT A ZERO''/&
        &'' STEP LENGTH'')',idg)
      CALL IVOUT(0,idum,'('' WHEN COL. NO. LEAVING=COL. NO. ENTERING, THE ENTERING &
        &VARIABLE MOVED''/&
        &'' TO ITS BOUND.  IT REMAINS NON-BASIC.''/&
        &'' WHEN COL. NO. OF BASIS EXCHANGED IS NEGATIVE, THE LEAVING''/&
        &'' VARIABLE IS AT ITS UPPER BOUND.'')',idg)
    ENDIF
  ENDIF
  !++  CODE FOR OUTPUT=NO IS INACTIVE
  !++  END
  !
  !     INITIALIZATION. SCALE DATA, NORMALIZE BOUNDS, FORM COLUMN
  !     CHECK SUMS, AND FORM INITIAL BASIS MATRIX.
  CALL DPINIT(Mrelas,Nvars,Costs,Bl,Bu,Ind,Primal,Info,Amat,Csc,costsc,&
    Colnrm,xlamda,anorm,Rhs,rhsnrm,Ibasis,Ibb,Imat,lopt)
  IF ( Info<0 ) GOTO 4600
  !
  nredc = 0
  npr004 = 200
  GOTO 2700
  200 CONTINUE
  IF ( .NOT.(singlr) ) THEN
    npr005 = 300
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PROCEDURE (COMPUTE ERROR IN DUAL AND PRIMAL SYSTEMS)
    ntries = 1
    GOTO 3000
  ELSE
    nerr = 23
    CALL XERMSG('SLATEC','DPLPMN',&
      'IN DSPLP,  A SINGULAR INITIAL BASIS WAS ENCOUNTERED.',nerr,&
      iopt)
    Info = -nerr
    GOTO 4600
  ENDIF
  300  npr006 = 400
  GOTO 4000
  400  npr007 = 500
  GOTO 2800
  500 CONTINUE
  IF ( .NOT.(usrbas) ) GOTO 700
  npr008 = 600
  GOTO 3100
  600 CONTINUE
  IF ( .NOT.feas ) THEN
    nerr = 24
    CALL XERMSG('SLATEC','DPLPMN',&
      'IN DSPLP, AN INFEASIBLE INITIAL BASIS WAS ENCOUNTERED.',&
      nerr,iopt)
    Info = -nerr
    GOTO 4600
  ENDIF
  700  itlp = 0
  !
  !     LAMDA HAS BEEN SET TO A CONSTANT, PERFORM PENALTY METHOD.
  npr009 = 800
  !     PROCEDURE (PERFORM SIMPLEX STEPS)
  npr013 = 2000
  GOTO 4100
  800  npr010 = 900
  GOTO 1900
  900  npr006 = 1000
  GOTO 4000
  1000 npr008 = 1100
  GOTO 3100
  1100 CONTINUE
  IF ( feas ) THEN
    !     CHECK IF ANY BASIC VARIABLES ARE STILL CLASSIFIED AS
    !     INFEASIBLE.  IF ANY ARE, THEN THIS MAY NOT YET BE AN
    !     OPTIMAL POINT.  THEREFORE SET LAMDA TO ZERO AND TRY
    !     TO PERFORM MORE SIMPLEX STEPS.
    i = 1
    n20046 = Mrelas
    DO WHILE ( (n20046-i)>=0 )
      IF ( Primal(i+Nvars)/=zero ) THEN
        xlamda = zero
        npr009 = 1700
        npr013 = 2000
        GOTO 4100
      ELSE
        i = i + 1
      ENDIF
    ENDDO
    GOTO 1700
  ELSE
    !
    !     SET LAMDA TO INFINITY BY SETTING COSTSC TO ZERO (SAVE THE VALUE OF
    !     COSTSC) AND PERFORM STANDARD PHASE-1.
    IF ( kprint>=2 ) CALL IVOUT(0,idum,'('' ENTER STANDARD PHASE-1'')',idg)
    scosts = costsc
    costsc = zero
    npr007 = 1200
    GOTO 2800
  ENDIF
  1200 npr009 = 1300
  npr013 = 2000
  GOTO 4100
  1300 npr010 = 1400
  GOTO 1900
  1400 npr006 = 1500
  GOTO 4000
  1500 npr008 = 1600
  GOTO 3100
  1600 CONTINUE
  IF ( feas ) THEN
    !
    !     SET LAMDA TO ZERO, COSTSC=SCOSTS, PERFORM STANDARD PHASE-2.
    IF ( kprint>1 ) CALL IVOUT(0,idum,'('' ENTER STANDARD PHASE-2'')',idg)
    xlamda = zero
    costsc = scosts
    npr009 = 1700
    npr013 = 2000
    GOTO 4100
  ENDIF
  !
  1700 npr011 = 1800
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE(RESCALE AND REARRANGE VARIABLES)
  !
  !     RESCALE THE DUAL VARIABLES.
  npr013 = 4300
  GOTO 4100
  1800 CONTINUE
  IF ( feas.AND.(.NOT.unbnd) ) THEN
    Info = 1
  ELSEIF ( (.NOT.feas).AND.(.NOT.unbnd) ) THEN
    nerr = 1
    CALL XERMSG('SLATEC','DPLPMN',&
      'IN DSPLP, THE PROBLEM APPEARS TO BE INFEASIBLE',nerr,iopt)
    Info = -nerr
  ELSEIF ( feas.AND.unbnd ) THEN
    nerr = 2
    CALL XERMSG('SLATEC','DPLPMN',&
      'IN DSPLP, THE PROBLEM APPEARS TO HAVE NO FINITE SOLUTION.',&
      nerr,iopt)
    Info = -nerr
  ELSEIF ( (.NOT.feas).AND.unbnd ) THEN
    nerr = 3
    CALL XERMSG('SLATEC','DPLPMN',&
      'IN DSPLP, THE PROBLEM APPEARS TO BE INFEASIBLE AND TO '//&
      'HAVE NO FINITE SOLN.',nerr,iopt)
    Info = -nerr
  ENDIF
  !
  IF ( Info==(-1).OR.Info==(-3) ) THEN
    size = DASUM(Nvars,Primal,1)*anorm
    size = size/DASUM(Nvars,Csc,1)
    size = size + DASUM(Mrelas,Primal(Nvars+1),1)
    i = 1
    n20058 = Nvars + Mrelas
    DO WHILE ( (n20058-i)>=0 )
      nx0066 = Ind(i)
      IF ( nx0066>=1.AND.nx0066<=4 ) THEN
        SELECT CASE (nx0066)
          CASE (2)
            IF ( size+ABS(Primal(i)-Bu(i))*factor/=size ) THEN
              IF ( Primal(i)>=Bu(i) ) Ind(i) = -4
            ENDIF
          CASE (3)
            IF ( size+ABS(Primal(i)-Bl(i))*factor/=size ) THEN
              IF ( Primal(i)<Bl(i) ) THEN
                Ind(i) = -4
              ELSEIF ( size+ABS(Primal(i)-Bu(i))*factor/=size ) THEN
                IF ( Primal(i)>Bu(i) ) Ind(i) = -4
              ENDIF
            ENDIF
          CASE (4)
          CASE DEFAULT
            IF ( size+ABS(Primal(i)-Bl(i))*factor/=size ) THEN
              IF ( Primal(i)<=Bl(i) ) Ind(i) = -4
            ENDIF
        END SELECT
      ENDIF
      i = i + 1
    ENDDO
  ENDIF
  !
  IF ( Info==(-2).OR.Info==(-3) ) THEN
    j = 1
    n20080 = Nvars
    DO WHILE ( (n20080-j)>=0 )
      IF ( Ibb(j)==0 ) THEN
        nx0091 = Ind(j)
        IF ( nx0091>=1.AND.nx0091<=4 ) THEN
          SELECT CASE (nx0091)
            CASE (2)
              Bl(j) = Bu(j)
              Ind(j) = -3
            CASE (3)
            CASE (4)
              Bl(j) = zero
              Bu(j) = zero
              Ind(j) = -3
            CASE DEFAULT
              Bu(j) = Bl(j)
              Ind(j) = -3
          END SELECT
        ENDIF
      ENDIF
      j = j + 1
    ENDDO
  ENDIF
  !++  CODE FOR OUTPUT=YES IS ACTIVE
  IF ( kprint<1 ) GOTO 4600
  npr012 = 4600
  !++  CODE FOR OUTPUT=NO IS INACTIVE
  !++  END
  GOTO 4500
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (COMPUTE RIGHT HAND SIDE)
  1900 Rhs(1) = zero
  CALL DCOPY(Mrelas,Rhs,0,Rhs,1)
  j = 1
  n20098 = Nvars + Mrelas
  DO WHILE ( (n20098-j)>=0 )
    nx0106 = Ind(j)
    IF ( nx0106>=1.AND.nx0106<=4 ) THEN
      SELECT CASE (nx0106)
        CASE (2)
          scalr = -Bu(j)
        CASE (3)
          scalr = -Bl(j)
        CASE (4)
          scalr = zero
        CASE DEFAULT
          scalr = -Bl(j)
      END SELECT
    ENDIF
    IF ( scalr==zero ) THEN
      j = j + 1
    ELSEIF ( j>Nvars ) THEN
      Rhs(j-Nvars) = Rhs(j-Nvars) - scalr
      j = j + 1
    ELSE
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,j)
        IF ( i>0 ) THEN
          Rhs(i) = Rhs(i) + aij*scalr
        ELSE
          j = j + 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  j = 1
  n20119 = Nvars + Mrelas
  DO WHILE ( (n20119-j)>=0 )
    scalr = zero
    IF ( Ind(j)==3.AND.MOD(Ibb(j),2)==0 ) scalr = Bu(j) - Bl(j)
    IF ( scalr==zero ) THEN
      j = j + 1
    ELSEIF ( j>Nvars ) THEN
      Rhs(j-Nvars) = Rhs(j-Nvars) + scalr
      j = j + 1
    ELSE
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,j)
        IF ( i>0 ) THEN
          Rhs(i) = Rhs(i) - aij*scalr
        ELSE
          j = j + 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  SELECT CASE(npr010)
    CASE(900)
      GOTO 900
    CASE(1400)
      GOTO 1400
  END SELECT
  2000 npr014 = 2100
  GOTO 3200
  2100 CONTINUE
  IF ( kprint>2 ) THEN
    CALL DVOUT(Mrelas,Duals,'('' BASIC (INTERNAL) DUAL SOLN.'')',idg)
    CALL DVOUT(Nvars+Mrelas,Rz,'('' REDUCED COSTS'')',idg)
  ENDIF
  npr015 = 2200
  GOTO 4200
  2200 CONTINUE
  IF ( .NOT.found ) THEN
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PROCEDURE (REDECOMPOSE BASIS MATRIX AND TRY AGAIN)
    IF ( redbas ) GOTO 3900
    npr004 = 3500
    GOTO 2700
  ENDIF
  2300 CONTINUE
  IF ( .NOT.(found) ) THEN
    SELECT CASE(npr009)
      CASE(800)
        GOTO 800
      CASE(1300)
        GOTO 1300
      CASE(1700)
        GOTO 1700
    END SELECT
  ELSE
    IF ( kprint>=3 ) CALL DVOUT(Mrelas,Ww,'('' SEARCH DIRECTION'')',idg)
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PROCEDURE (CHOOSE VARIABLE TO LEAVE BASIS)
    CALL DPLPFL(Mrelas,Nvars,ienter,ileave,Ibasis,Ind,Ibb,theta,dirnrm,&
      rprnrm,Csc,Ww,Bl,Bu,Erp,Rprim,Primal,finite,zerolv)
    IF ( .NOT.(finite) ) THEN
      unbnd = .TRUE.
      Ibb(Ibasis(ienter)) = 0
    ELSE
      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !     PROCEDURE (MAKE MOVE AND UPDATE)
      CALL DPLPMU(Mrelas,Nvars,Lmx,Lbm,nredc,Info,ienter,ileave,iopt,npp,&
        jstrt,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,anorm,eps,uu,gg,&
        rprnrm,erdnrm,dulnrm,theta,costsc,xlamda,rhsnrm,Amat,&
        Basmat,Csc,Wr,Rprim,Ww,Bu,Bl,Rhs,Erd,Erp,Rz,Rg,Colnrm,&
        Costs,Primal,Duals,singlr,redbas,zerolv,stpedg)
      IF ( Info==(-26) ) GOTO 4600
      !++  CODE FOR OUTPUT=YES IS ACTIVE
      IF ( kprint>=2 ) THEN
        ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !     PROCEDURE (PRINT ITERATION SUMMARY)
        idum(1) = itlp + 1
        CALL IVOUT(1,idum,'(''0ITERATION NUMBER'')',idg)
        idum(1) = Ibasis(ABS(ileave))
        CALL IVOUT(1,idum,'('' INDEX OF VARIABLE ENTERING THE BASIS'')',idg)
        idum(1) = ileave
        CALL IVOUT(1,idum,'('' COLUMN OF THE BASIS EXCHANGED'')',idg)
        idum(1) = Ibasis(ienter)
        CALL IVOUT(1,idum,'('' INDEX OF VARIABLE LEAVING THE BASIS'')',idg)
        rdum(1) = theta
        CALL DVOUT(1,rdum,'('' LENGTH OF THE EXCHANGE STEP'')',idg)
        IF ( kprint>=3 ) THEN
          CALL DVOUT(Mrelas,Rprim,'('' BASIC (INTERNAL) PRIMAL SOLN.'')',&
            idg)
          CALL IVOUT(Nvars+Mrelas,Ibasis,&
            '('' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE BASIC.'')'&
            ,idg)
          CALL IVOUT(Nvars+Mrelas,Ibb,'('' IBB ARRAY'')',idg)
          CALL DVOUT(Mrelas,Rhs,'('' TRANSLATED RHS'')',idg)
          CALL DVOUT(Mrelas,Duals,'('' BASIC (INTERNAL) DUAL SOLN.'')',idg)
        ENDIF
        !++  CODE FOR OUTPUT=NO IS INACTIVE
        !++  END
      ENDIF
      npr005 = 2400
      ntries = 1
      GOTO 3000
    ENDIF
  ENDIF
  2400 itlp = itlp + 1
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (CHECK AND RETURN WITH EXCESS ITERATIONS)
  IF ( itlp<=mxitlp ) THEN
    npr015 = 2200
    GOTO 4200
  ELSE
    nerr = 25
    npr011 = 3300
    npr013 = 4300
    GOTO 4100
  ENDIF
  2500 lpr1 = lpr + 1
  READ (isave) (Amat(i),i=lpr1,Lmx), (Imat(i),i=lpr1,Lmx)
  CALL DPRWPG(key,ipage,lpg,Amat,Imat)
  np = Imat(Lmx-1)
  ipage = ipage + 1
  IF ( np>=0 ) GOTO 2500
  nparm = Nvars + Mrelas
  READ (isave) (Ibasis(i),i=1,nparm)
  REWIND isave
  GOTO 100
  2600 CALL DPRWPG(key,ipage,lpg,Amat,Imat)
  lpr1 = lpr + 1
  WRITE (isave) (Amat(i),i=lpr1,Lmx), (Imat(i),i=lpr1,Lmx)
  np = Imat(Lmx-1)
  ipage = ipage + 1
  IF ( np>=0 ) GOTO 2600
  nparm = Nvars + Mrelas
  WRITE (isave) (Ibasis(i),i=1,nparm)
  ENDFILE isave
  IF ( Imat(Lmx-1)/=(-1) ) CALL SCLOSM(ipagef)
  GOTO 99999
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (DECOMPOSE BASIS MATRIX)
  !++  CODE FOR OUTPUT=YES IS ACTIVE
  2700 CONTINUE
  IF ( kprint>=2 ) CALL IVOUT(Mrelas,Ibasis,&
    '('' SUBSCRIPTS OF BASIC VARIABLES DURING REDECOMPOSITION'')'&
    ,idg)
  !++  CODE FOR OUTPUT=NO IS INACTIVE
  !++  END
  !
  !     SET RELATIVE PIVOTING FACTOR FOR USE IN LA05 () PACKAGE.
  uu = 0.1
  CALL DPLPDM(Mrelas,Nvars,Lmx,Lbm,nredc,Info,iopt,Ibasis,Imat,Ibrc,Ipr,Iwr,&
    Ind,Ibb,anorm,eps,uu,gg,Amat,Basmat,Csc,Wr,singlr,redbas)
  IF ( Info<0 ) GOTO 4600
  SELECT CASE(npr004)
    CASE(200)
      GOTO 200
    CASE(2900)
      GOTO 2900
    CASE(3500)
      GOTO 3500
  END SELECT
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (CLASSIFY VARIABLES)
  !
  !     DEFINE THE CLASSIFICATION OF THE BASIC VARIABLES
  !     -1 VIOLATES LOWER BOUND, 0 FEASIBLE, +1 VIOLATES UPPER BOUND.
  !     (THIS INFO IS STORED IN PRIMAL(NVARS+1)-PRIMAL(NVARS+MRELAS))
  !     TRANSLATE VARIABLE TO ITS UPPER BOUND, IF .GT. UPPER BOUND
  2800 Primal(Nvars+1) = zero
  CALL DCOPY(Mrelas,Primal(Nvars+1),0,Primal(Nvars+1),1)
  i = 1
  n20172 = Mrelas
  DO WHILE ( (n20172-i)>=0 )
    j = Ibasis(i)
    IF ( Ind(j)/=4 ) THEN
      IF ( Rprim(i)<zero ) THEN
        Primal(i+Nvars) = -one
      ELSEIF ( Ind(j)==3 ) THEN
        upbnd = Bu(j) - Bl(j)
        IF ( j<=Nvars ) upbnd = upbnd/Csc(j)
        IF ( Rprim(i)>upbnd ) THEN
          Rprim(i) = Rprim(i) - upbnd
          IF ( j>Nvars ) THEN
            Rhs(j-Nvars) = Rhs(j-Nvars) + upbnd
          ELSE
            k = 0
            DO
              CALL DPNNZR(k,aij,iplace,Amat,Imat,j)
              IF ( k<=0 ) EXIT
              Rhs(k) = Rhs(k) - upbnd*aij*Csc(j)
            ENDDO
          ENDIF
          Primal(i+Nvars) = one
        ENDIF
      ENDIF
    ENDIF
    i = i + 1
  ENDDO
  SELECT CASE(npr007)
    CASE(500)
      GOTO 500
    CASE(1200)
      GOTO 1200
  END SELECT
  2900 ntries = ntries + 1
  3000 CONTINUE
  IF ( (2-ntries)>=0 ) THEN
    CALL DPLPCE(Mrelas,Nvars,Lmx,Lbm,itlp,itbrc,Ibasis,Imat,Ibrc,Ipr,Iwr,&
      Ind,Ibb,erdnrm,eps,tune,gg,Amat,Basmat,Csc,Wr,Ww,Primal,Erd,&
      Erp,singlr,redbas)
    IF ( .NOT.singlr ) THEN
      !++  CODE FOR OUTPUT=YES IS ACTIVE
      IF ( kprint>=3 ) THEN
        CALL DVOUT(Mrelas,Erp,'('' EST. ERROR IN PRIMAL COMPS.'')',idg)
        !++  CODE FOR OUTPUT=NO IS INACTIVE
        !++  END
        CALL DVOUT(Mrelas,Erd,'('' EST. ERROR IN DUAL COMPS.'')',idg)
      ENDIF
      SELECT CASE(npr005)
        CASE(300)
          GOTO 300
        CASE(2400)
          GOTO 2400
        CASE(3600)
          GOTO 3600
      END SELECT
    ELSEIF ( ntries/=2 ) THEN
      npr004 = 2900
      GOTO 2700
    ENDIF
  ENDIF
  nerr = 26
  CALL XERMSG('SLATEC','DPLPMN',&
    'IN DSPLP, MOVED TO A SINGULAR POINT. THIS SHOULD NOT HAPPEN.'&
    ,nerr,iopt)
  Info = -nerr
  GOTO 4600
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (CHECK FEASIBILITY)
  !
  !     SEE IF NEARBY FEASIBLE POINT SATISFIES THE CONSTRAINT
  !     EQUATIONS.
  !
  !     COPY RHS INTO WW(*), THEN UPDATE WW(*).
  3100 CALL DCOPY(Mrelas,Rhs,1,Ww,1)
  j = 1
  n20206 = Mrelas
  DO WHILE ( (n20206-j)>=0 )
    ibas = Ibasis(j)
    xval = Rprim(j)
    !
    !     ALL VARIABLES BOUNDED BELOW HAVE ZERO AS THAT BOUND.
    IF ( Ind(ibas)<=3 ) xval = MAX(zero,xval)
    !
    !     IF THE VARIABLE HAS AN UPPER BOUND, COMPUTE THAT BOUND.
    IF ( Ind(ibas)==3 ) THEN
      upbnd = Bu(ibas) - Bl(ibas)
      IF ( ibas<=Nvars ) upbnd = upbnd/Csc(ibas)
      xval = MIN(upbnd,xval)
    ENDIF
    !
    !     SUBTRACT XVAL TIMES COLUMN VECTOR FROM RIGHT-HAND SIDE IN WW(*)
    IF ( xval==zero ) THEN
      j = j + 1
    ELSEIF ( ibas>Nvars ) THEN
      IF ( Ind(ibas)/=2 ) THEN
        Ww(ibas-Nvars) = Ww(ibas-Nvars) + xval
      ELSE
        Ww(ibas-Nvars) = Ww(ibas-Nvars) - xval
      ENDIF
      j = j + 1
    ELSE
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,ibas)
        IF ( i>0 ) THEN
          Ww(i) = Ww(i) - xval*aij*Csc(ibas)
        ELSE
          j = j + 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
  !   COMPUTE NORM OF DIFFERENCE AND CHECK FOR FEASIBILITY.
  resnrm = DASUM(Mrelas,Ww,1)
  feas = resnrm<=tolls*(rprnrm*anorm+rhsnrm)
  !
  !     TRY AN ABSOLUTE ERROR TEST IF THE RELATIVE TEST FAILS.
  IF ( .NOT.feas ) feas = resnrm<=tolabs
  IF ( feas ) THEN
    Primal(Nvars+1) = zero
    CALL DCOPY(Mrelas,Primal(Nvars+1),0,Primal(Nvars+1),1)
  ENDIF
  SELECT CASE(npr008)
    CASE(600)
      GOTO 600
    CASE(1100)
      GOTO 1100
    CASE(1600)
      GOTO 1600
  END SELECT
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (INITIALIZE REDUCED COSTS AND STEEPEST EDGE WEIGHTS)
  3200 CALL DPINCW(Mrelas,Nvars,Lmx,Lbm,npp,jstrt,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,&
    Ibb,costsc,gg,erdnrm,dulnrm,Amat,Basmat,Csc,Wr,Ww,Rz,Rg,Costs,&
    Colnrm,Duals,stpedg)
  !
  SELECT CASE(npr014)
    CASE(2100)
      GOTO 2100
    CASE(3900)
      GOTO 3900
  END SELECT
  !++  CODE FOR OUTPUT=YES IS ACTIVE
  3300 CONTINUE
  IF ( kprint>=1 ) THEN
    npr012 = 3400
    GOTO 4500
  ENDIF
  !++  CODE FOR OUTPUT=NO IS INACTIVE
  !++  END
  3400 idum(1) = 0
  IF ( savedt ) idum(1) = isave
  WRITE (xern1,'(I8)') mxitlp
  WRITE (xern2,'(I8)') idum(1)
  CALL XERMSG('SLATEC','DPLPMN','IN DSPLP, MAX ITERATIONS = '//xern1//&
    ' TAKEN.  UP-TO-DATE RESULTS SAVED ON FILE NO. '//xern2//&
    '.   IF FILE NO. = 0, NO SAVE.',nerr,iopt)
  Info = -nerr
  GOTO 4600
  3500 npr005 = 3600
  ntries = 1
  GOTO 3000
  3600 npr006 = 3700
  GOTO 4000
  3700 npr013 = 3800
  GOTO 4100
  3800 npr014 = 3900
  GOTO 3200
  !
  !     ERASE NON-CYCLING MARKERS NEAR COMPLETION.
  3900 i = Mrelas + 1
  n20247 = Mrelas + Nvars
  DO WHILE ( (n20247-i)>=0 )
    Ibasis(i) = ABS(Ibasis(i))
    i = i + 1
  ENDDO
  npr015 = 2300
  GOTO 4200
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (COMPUTE NEW PRIMAL)
  !
  !     COPY RHS INTO WW(*), SOLVE SYSTEM.
  4000 CALL DCOPY(Mrelas,Rhs,1,Ww,1)
  trans = .FALSE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,gg,Ww,trans)
  CALL DCOPY(Mrelas,Ww,1,Rprim,1)
  rprnrm = DASUM(Mrelas,Rprim,1)
  SELECT CASE(npr006)
    CASE(400)
      GOTO 400
    CASE(1000)
      GOTO 1000
    CASE(1500)
      GOTO 1500
    CASE(3700)
      GOTO 3700
    CASE(4400)
      GOTO 4400
  END SELECT
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (COMPUTE NEW DUALS)
  !
  !     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*).
  4100 i = 1
  n20252 = Mrelas
  DO WHILE ( (n20252-i)>=0 )
    j = Ibasis(i)
    IF ( j>Nvars ) THEN
      Duals(i) = xlamda*Primal(i+Nvars)
    ELSE
      Duals(i) = costsc*Costs(j)*Csc(j) + xlamda*Primal(i+Nvars)
    ENDIF
    i = i + 1
  ENDDO
  !
  trans = .TRUE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,gg,Duals,trans)
  dulnrm = DASUM(Mrelas,Duals,1)
  SELECT CASE(npr013)
    CASE(2000)
      GOTO 2000
    CASE(3800)
      GOTO 3800
    CASE(4300)
      GOTO 4300
  END SELECT
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (FIND VARIABLE TO ENTER BASIS AND GET SEARCH DIRECTION)
  4200 CALL DPLPFE(Mrelas,Nvars,Lmx,Lbm,ienter,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,&
    erdnrm,eps,gg,dulnrm,dirnrm,Amat,Basmat,Csc,Wr,Ww,Bl,Bu,Rz,Rg,&
    Colnrm,Duals,found)
  SELECT CASE(npr015)
    CASE(2200)
      GOTO 2200
    CASE(2300)
      GOTO 2300
  END SELECT
  4300 CONTINUE
  IF ( costsc==zero ) THEN
    npr006 = 4400
    GOTO 4000
  ELSE
    i = 1
    n20271 = Mrelas
    DO WHILE ( (n20271-i)>=0 )
      Duals(i) = Duals(i)/costsc
      i = i + 1
    ENDDO
    npr006 = 4400
    GOTO 4000
  ENDIF
  !
  !     REAPPLY COLUMN SCALING TO PRIMAL.
  4400 i = 1
  n20276 = Mrelas
  DO WHILE ( (n20276-i)>=0 )
    j = Ibasis(i)
    IF ( j<=Nvars ) THEN
      scalr = Csc(j)
      IF ( Ind(j)==2 ) scalr = -scalr
      Rprim(i) = Rprim(i)*scalr
    ENDIF
    i = i + 1
  ENDDO
  !
  !     REPLACE TRANSLATED BASIC VARIABLES INTO ARRAY PRIMAL(*)
  Primal(1) = zero
  CALL DCOPY(Nvars+Mrelas,Primal,0,Primal,1)
  j = 1
  n20283 = Nvars + Mrelas
  DO WHILE ( (n20283-j)>=0 )
    ibas = ABS(Ibasis(j))
    xval = zero
    IF ( j<=Mrelas ) xval = Rprim(j)
    IF ( Ind(ibas)==1 ) xval = xval + Bl(ibas)
    IF ( Ind(ibas)==2 ) xval = Bu(ibas) - xval
    IF ( Ind(ibas)==3 ) THEN
      IF ( MOD(Ibb(ibas),2)==0 ) xval = Bu(ibas) - Bl(ibas) - xval
      xval = xval + Bl(ibas)
    ENDIF
    Primal(ibas) = xval
    j = j + 1
  ENDDO
  !
  !     COMPUTE DUALS FOR INDEPENDENT VARIABLES WITH BOUNDS.
  !     OTHER ENTRIES ARE ZERO.
  j = 1
  n20290 = Nvars
  DO WHILE ( (n20290-j)>=0 )
    rzj = zero
    IF ( Ibb(j)>zero.AND.Ind(j)/=4 ) THEN
      rzj = Costs(j)
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,j)
        IF ( i<=0 ) EXIT
        rzj = rzj - aij*Duals(i)
      ENDDO
    ENDIF
    Duals(Mrelas+j) = rzj
    j = j + 1
  ENDDO
  SELECT CASE(npr011)
    CASE(1800)
      GOTO 1800
    CASE(3300)
      GOTO 3300
  END SELECT
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (PRINT SUMMARY)
  4500 idum(1) = Info
  CALL IVOUT(1,idum,'('' THE OUTPUT VALUE OF INFO IS'')',idg)
  IF ( .NOT.(minprb) ) THEN
    CALL IVOUT(0,idum,'('' THIS IS A MAXIMIZATION PROBLEM.'')',idg)
  ELSE
    CALL IVOUT(0,idum,'('' THIS IS A MINIMIZATION PROBLEM.'')',idg)
  ENDIF
  IF ( .NOT.(stpedg) ) THEN
    CALL IVOUT(0,idum,'('' MINIMUM REDUCED COST PRICING WAS USED.'')',idg)
  ELSE
    CALL IVOUT(0,idum,'('' STEEPEST EDGE PRICING WAS USED.'')',idg)
  ENDIF
  rdum(1) = DDOT(Nvars,Costs,1,Primal,1)
  CALL DVOUT(1,rdum,'('' OUTPUT VALUE OF THE OBJECTIVE FUNCTION'')',idg)
  CALL DVOUT(Nvars+Mrelas,Primal,&
    '('' THE OUTPUT INDEPENDENT AND DEPENDENT VARIABLES'')',idg)
  CALL DVOUT(Mrelas+Nvars,Duals,'('' THE OUTPUT DUAL VARIABLES'')',idg)
  CALL IVOUT(Nvars+Mrelas,Ibasis,&
    '('' VARIABLE INDICES IN POSITIONS 1-MRELAS ARE BASIC.'')',idg)
  idum(1) = itlp
  CALL IVOUT(1,idum,'('' NO. OF ITERATIONS'')',idg)
  idum(1) = nredc
  CALL IVOUT(1,idum,'('' NO. OF FULL REDECOMPS'')',idg)
  SELECT CASE(npr012)
    CASE(3400)
      GOTO 3400
    CASE(4600)
      GOTO 4600
  END SELECT
  !++  CODE FOR OUTPUT=NO IS INACTIVE
  !++  END
  ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !     PROCEDURE (RETURN TO USER)
  4600 CONTINUE
  IF ( .NOT.(savedt) ) THEN
    IF ( Imat(Lmx-1)/=(-1) ) CALL SCLOSM(ipagef)
  ELSE
    ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !     PROCEDURE (SAVE DATA ON FILE ISAVE)
    !
    !     SOME PAGES MAY NOT BE WRITTEN YET.
    IF ( Amat(Lmx)==one ) THEN
      Amat(Lmx) = zero
      key = 2
      ipage = ABS(Imat(Lmx-1))
      CALL DPRWPG(key,ipage,lpg,Amat,Imat)
    ENDIF
    !
    !     FORCE PAGE FILE TO BE OPENED ON RESTARTS.
    key = Amat(4)
    Amat(4) = zero
    lpr = Nvars + 4
    WRITE (isave) (Amat(i),i=1,lpr), (Imat(i),i=1,lpr)
    Amat(4) = key
    ipage = 1
    key = 1
    GOTO 2600
  ENDIF
  !
  !     THIS TEST IS THERE ONLY TO AVOID DIAGNOSTICS ON SOME FORTRAN
  !     COMPILERS.
  99999 CONTINUE
  END SUBROUTINE DPLPMN
