!*==DPLPMU.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DPLPMU
SUBROUTINE DPLPMU(Mrelas,Nvars,Lmx,Lbm,Nredc,Info,Ienter,Ileave,Iopt,Npp,&
    Jstrt,Ibasis,Imat,Ibrc,Ipr,Iwr,Ind,Ibb,Anorm,Eps,Uu,Gg,&
    Rprnrm,Erdnrm,Dulnrm,Theta,Costsc,Xlamda,Rhsnrm,Amat,&
    Basmat,Csc,Wr,Rprim,Ww,Bu,Bl,Rhs,Erd,Erp,Rz,Rg,Colnrm,&
    Costs,Primal,Duals,Singlr,Redbas,Zerolv,Stpedg)
  IMPLICIT NONE
  !*--DPLPMU9
  !*** Start of declarations inserted by SPAG
  INTEGER i , ibas , IDLOC , Ienter , ihi , il1 , Ileave , ilow , Info , &
    Iopt , ipage , iplace , iu1 , j , Jstrt , k , key , Lbm , Lmx , &
    lpg
  INTEGER Mrelas , n20002 , n20018 , n20121 , nerr , nnegrc , Npp , npr001 , &
    npr003 , Nredc , Nvars
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DPLPMU
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SPLPMU-S, DPLPMU-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/,
  !     /SASUM/DASUM/,/SCOPY/DCOPY/,/SDOT/DDOT/,
  !     /.E0/.D0/
  !
  !     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE
  !     TASKS OF UPDATING THE PRIMAL SOLUTION, EDGE WEIGHTS, REDUCED
  !     COSTS, AND MATRIX DECOMPOSITION.
  !     IT IS THE MAIN PART OF THE PROCEDURE (MAKE MOVE AND UPDATE).
  !
  !     REVISED 821122-1100
  !     REVISED YYMMDD
  !
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  DASUM, DCOPY, DDOT, DPLPDM, DPNNZR, DPRWPG, IDLOC,
  !                    LA05BD, LA05CD, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   890606  Removed unused COMMON block LA05DD.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DPLPMU
  INTEGER Ibasis(*) , Imat(*) , Ibrc(Lbm,2) , Ipr(*) , Iwr(*) , Ind(*) , &
    Ibb(*)
  DOUBLE PRECISION aij , alpha , Anorm , Costsc , Erdnrm , Dulnrm , Eps , &
    gamma , Gg , gq , one , Rprnrm , rzj , scalr , Theta , &
    two , Uu , wp , Xlamda , Rhsnrm , zero , Amat(*) , &
    Basmat(*) , Csc(*) , Wr(*) , Rprim(*) , Ww(*) , Bu(*) , &
    Bl(*) , Rhs(*) , Erd(*) , Erp(*) , Rz(*) , Rg(*) , &
    Costs(*) , Primal(*) , Duals(*) , Colnrm(*) , rcost , &
    DASUM , DDOT , cnorm
  LOGICAL Singlr , Redbas , pagepl , trans , Zerolv , Stpedg
  !
  !***FIRST EXECUTABLE STATEMENT  DPLPMU
  zero = 0.D0
  one = 1.D0
  two = 2.D0
  lpg = Lmx - (Nvars+4)
  !
  !     UPDATE THE PRIMAL SOLUTION WITH A MULTIPLE OF THE SEARCH
  !     DIRECTION.
  i = 1
  n20002 = Mrelas
  DO WHILE ( (n20002-i)>=0 )
    Rprim(i) = Rprim(i) - Theta*Ww(i)
    i = i + 1
  ENDDO
  !
  !     IF EJECTED VARIABLE IS LEAVING AT AN UPPER BOUND,  THEN
  !     TRANSLATE RIGHT HAND SIDE.
  IF ( Ileave>=0 ) GOTO 200
  ibas = Ibasis(ABS(Ileave))
  scalr = Rprim(ABS(Ileave))
  npr001 = 100
  GOTO 1800
  100  Ibb(ibas) = ABS(Ibb(ibas)) + 1
  !
  !     IF ENTERING VARIABLE IS RESTRICTED TO ITS UPPER BOUND, TRANSLATE
  !     RIGHT HAND SIDE.  IF THE VARIABLE DECREASED FROM ITS UPPER
  !     BOUND, A SIGN CHANGE IS REQUIRED IN THE TRANSLATION.
  200  IF ( Ienter/=Ileave ) THEN
  ibas = Ibasis(Ienter)
  !
  !     IF ENTERING VARIABLE IS DECREASING FROM ITS UPPER BOUND,
  !     COMPLEMENT ITS PRIMAL VALUE.
  IF ( Ind(ibas)/=3.OR.MOD(Ibb(ibas),2)/=0 ) GOTO 500
  scalr = -(Bu(ibas)-Bl(ibas))
  IF ( ibas<=Nvars ) scalr = scalr/Csc(ibas)
  npr001 = 400
  GOTO 1800
ELSE
  ibas = Ibasis(Ienter)
  scalr = Theta
  IF ( MOD(Ibb(ibas),2)==0 ) scalr = -scalr
  npr001 = 300
  GOTO 1800
ENDIF
300  Ibb(ibas) = Ibb(ibas) + 1
GOTO 600
400  Theta = -scalr - Theta
Ibb(ibas) = Ibb(ibas) + 1
500  Rprim(ABS(Ileave)) = Theta
Ibb(ibas) = -ABS(Ibb(ibas))
i = Ibasis(ABS(Ileave))
Ibb(i) = ABS(Ibb(i))
IF ( Primal(ABS(Ileave)+Nvars)>zero ) Ibb(i) = Ibb(i) + 1
!
!     INTERCHANGE COLUMN POINTERS TO NOTE EXCHANGE OF COLUMNS.
600  ibas = Ibasis(Ienter)
Ibasis(Ienter) = Ibasis(ABS(Ileave))
Ibasis(ABS(Ileave)) = ibas
!
!     IF VARIABLE WAS EXCHANGED AT A ZERO LEVEL, MARK IT SO THAT
!     IT CAN'T BE BROUGHT BACK IN.  THIS IS TO HELP PREVENT CYCLING.
IF ( Zerolv ) Ibasis(Ienter) = -ABS(Ibasis(Ienter))
Rprnrm = MAX(Rprnrm,DASUM(Mrelas,Rprim,1))
k = 1
n20018 = Mrelas
700  DO WHILE ( (n20018-k)>=0 )
  !
  !     SEE IF VARIABLES THAT WERE CLASSIFIED AS INFEASIBLE HAVE NOW
  !     BECOME FEASIBLE.  THIS MAY REQUIRED TRANSLATING UPPER BOUNDED
  !     VARIABLES.
  IF ( Primal(k+Nvars)==zero.OR.ABS(Rprim(k))>Rprnrm*Erp(k) ) THEN
    k = k + 1
  ELSE
    IF ( Primal(k+Nvars)<=zero ) GOTO 900
    ibas = Ibasis(k)
    scalr = -(Bu(ibas)-Bl(ibas))
    IF ( ibas<=Nvars ) scalr = scalr/Csc(ibas)
    npr001 = 800
    GOTO 1800
  ENDIF
ENDDO
!
!     UPDATE REDUCED COSTS, EDGE WEIGHTS, AND MATRIX DECOMPOSITION.
IF ( Ienter==Ileave ) THEN
  !
  !     THIS IS NECESSARY ONLY FOR PRINTING OF INTERMEDIATE RESULTS.
  npr003 = 1700
  GOTO 1900
ELSE
  !
  !     THE INCOMING VARIABLE IS ALWAYS CLASSIFIED AS FEASIBLE.
  Primal(ABS(Ileave)+Nvars) = zero
  !
  wp = Ww(ABS(Ileave))
  gq = DDOT(Mrelas,Ww,1,Ww,1) + one
  !
  !     COMPUTE INVERSE (TRANSPOSE) TIMES SEARCH DIRECTION.
  trans = .TRUE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
  !
  !     UPDATE THE MATRIX DECOMPOSITION.  COL. ABS(ILEAVE) IS LEAVING.
  !     THE ARRAY DUALS(*) CONTAINS INTERMEDIATE RESULTS FOR THE
  !     INCOMING COLUMN.
  CALL LA05CD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Duals,Gg,Uu,ABS(Ileave))
  Redbas = .FALSE.
  IF ( Gg>=zero ) GOTO 1000
  !
  !     REDECOMPOSE BASIS MATRIX WHEN AN ERROR RETURN FROM
  !     LA05CD( ) IS NOTED.  THIS WILL PROBABLY BE DUE TO
  !     SPACE BEING EXHAUSTED, GG=-7.
  CALL DPLPDM(Mrelas,Nvars,Lmx,Lbm,Nredc,Info,Iopt,Ibasis,Imat,Ibrc,Ipr,&
    Iwr,Ind,Ibb,Anorm,Eps,Uu,Gg,Amat,Basmat,Csc,Wr,Singlr,&
    Redbas)
  IF ( .NOT.(Singlr) ) THEN
    !     PROCEDURE (COMPUTE NEW PRIMAL)
    !
    !     COPY RHS INTO WW(*), SOLVE SYSTEM.
    CALL DCOPY(Mrelas,Rhs,1,Ww,1)
    trans = .FALSE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
    CALL DCOPY(Mrelas,Ww,1,Rprim,1)
    Rprnrm = DASUM(Mrelas,Rprim,1)
    GOTO 1000
  ELSE
    nerr = 26
    CALL XERMSG('SLATEC','DPLPMU',&
      'IN DSPLP, MOVED TO A SINGULAR POINT. THIS SHOULD NOT HAPPEN.'&
      ,nerr,Iopt)
    Info = -nerr
    RETURN
  ENDIF
ENDIF
800  Rprim(k) = -scalr
Rprnrm = Rprnrm - scalr
900  Primal(k+Nvars) = zero
k = k + 1
GOTO 700
!
!     IF STEEPEST EDGE PRICING IS USED, UPDATE REDUCED COSTS
!     AND EDGE WEIGHTS.
1000 IF ( .NOT.(Stpedg) ) THEN
!
!     COMPUTE THE UPDATED DUALS IN DUALS(*).
npr003 = 1300
ELSE
!
!     COMPUTE COL. ABS(ILEAVE) OF THE NEW INVERSE (TRANSPOSE) MATRIX
!     HERE ABS(ILEAVE) POINTS TO THE EJECTED COLUMN.
!     USE ERD(*) FOR TEMP. STORAGE.
CALL DCOPY(Mrelas,zero,0,Erd,1)
Erd(ABS(Ileave)) = one
trans = .TRUE.
CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Erd,trans)
!
!     COMPUTE UPDATED DUAL VARIABLES IN DUALS(*).
npr003 = 1100
ENDIF
GOTO 1900
!
!     COMPUTE THE DOT PRODUCT OF COL. J OF THE NEW INVERSE (TRANSPOSE)
!     WITH EACH NON-BASIC COLUMN.  ALSO COMPUTE THE DOT PRODUCT OF THE
!     INVERSE (TRANSPOSE) OF NON-UPDATED MATRIX (TIMES) THE
!     SEARCH DIRECTION WITH EACH NON-BASIC COLUMN.
!     RECOMPUTE REDUCED COSTS.
1100 pagepl = .TRUE.
CALL DCOPY(Nvars+Mrelas,zero,0,Rz,1)
nnegrc = 0
j = Jstrt
1200 IF ( Ibb(j)<=0 ) THEN
pagepl = .TRUE.
Rg(j) = one
!
!     NONBASIC INDEPENDENT VARIABLES (COLUMN IN SPARSE MATRIX STORAGE)
ELSEIF ( j>Nvars ) THEN
pagepl = .TRUE.
scalr = -one
IF ( Ind(j)==2 ) scalr = one
i = j - Nvars
alpha = scalr*Erd(i)
Rz(j) = -scalr*Duals(i)
gamma = scalr*Ww(i)
Rg(j) = MAX(Rg(j)-two*alpha*gamma+alpha**2*gq,one+alpha**2)
ELSE
rzj = Costs(j)*Costsc
alpha = zero
gamma = zero
!
!     COMPUTE THE DOT PRODUCT OF THE SPARSE MATRIX NONBASIC COLUMNS
!     WITH THREE VECTORS INVOLVED IN THE UPDATING STEP.
IF ( j/=1 ) THEN
ilow = Imat(j+3) + 1
ELSE
ilow = Nvars + 5
ENDIF
IF ( .NOT.(pagepl) ) THEN
il1 = ihi + 1
ELSE
il1 = IDLOC(ilow,Amat,Imat)
IF ( il1>=Lmx-1 ) THEN
  ilow = ilow + 2
  il1 = IDLOC(ilow,Amat,Imat)
ENDIF
ipage = ABS(Imat(Lmx-1))
ENDIF
ihi = Imat(j+4) - (ilow-il1)
DO
iu1 = MIN(Lmx-2,ihi)
IF ( il1>iu1 ) EXIT
DO i = il1 , iu1
  rzj = rzj - Amat(i)*Duals(Imat(i))
  alpha = alpha + Amat(i)*Erd(Imat(i))
  gamma = gamma + Amat(i)*Ww(Imat(i))
ENDDO
IF ( ihi<=Lmx-2 ) EXIT
ipage = ipage + 1
key = 1
CALL DPRWPG(key,ipage,lpg,Amat,Imat)
il1 = Nvars + 5
ihi = ihi - lpg
ENDDO
pagepl = ihi==(Lmx-2)
Rz(j) = rzj*Csc(j)
alpha = alpha*Csc(j)
gamma = gamma*Csc(j)
!
!     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
Rg(j) = MAX(Rg(j)-two*alpha*gamma+alpha**2*gq,one+alpha**2)
ENDIF
!
rcost = Rz(j)
IF ( MOD(Ibb(j),2)==0 ) rcost = -rcost
IF ( Ind(j)==3 ) THEN
IF ( Bu(j)==Bl(j) ) rcost = zero
ENDIF
IF ( Ind(j)==4 ) rcost = -ABS(rcost)
cnorm = one
IF ( j<=Nvars ) cnorm = Colnrm(j)
IF ( rcost+Erdnrm*Dulnrm*cnorm<zero ) nnegrc = nnegrc + 1
j = MOD(j,Mrelas+Nvars) + 1
IF ( nnegrc<Npp.AND.j/=Jstrt ) GOTO 1200
Jstrt = j
!
!     UPDATE THE EDGE WEIGHT FOR THE EJECTED VARIABLE.
Rg(ABS(Ibasis(Ienter))) = gq/wp**2
!
!     IF MINIMUM REDUCED COST (DANTZIG) PRICING IS USED,
!     CALCULATE THE NEW REDUCED COSTS.
GOTO 1700
1300 CALL DCOPY(Nvars+Mrelas,zero,0,Rz,1)
nnegrc = 0
j = Jstrt
pagepl = .TRUE.
!
1400 IF ( Ibb(j)<=0 ) THEN
pagepl = .TRUE.
GOTO 1600
!
!     NONBASIC INDEPENDENT VARIABLE (COLUMN IN SPARSE MATRIX STORAGE)
ELSEIF ( j>Nvars ) THEN
pagepl = .TRUE.
scalr = -one
IF ( Ind(j)==2 ) scalr = one
i = j - Nvars
Rz(j) = -scalr*Duals(i)
GOTO 1600
ELSE
Rz(j) = Costs(j)*Costsc
IF ( j/=1 ) THEN
ilow = Imat(j+3) + 1
ELSE
ilow = Nvars + 5
ENDIF
IF ( .NOT.(pagepl) ) THEN
il1 = ihi + 1
ELSE
il1 = IDLOC(ilow,Amat,Imat)
IF ( il1>=Lmx-1 ) THEN
ilow = ilow + 2
il1 = IDLOC(ilow,Amat,Imat)
ENDIF
ipage = ABS(Imat(Lmx-1))
ENDIF
ihi = Imat(j+4) - (ilow-il1)
ENDIF
1500 iu1 = MIN(Lmx-2,ihi)
IF ( iu1>=il1.AND.MOD(iu1-il1,2)==0 ) THEN
Rz(j) = Rz(j) - Amat(il1)*Duals(Imat(il1))
il1 = il1 + 1
ENDIF
IF ( il1<=iu1 ) THEN
!
!     UNROLL THE DOT PRODUCT LOOP TO A DEPTH OF TWO.  (THIS IS DONE
!     FOR INCREASED EFFICIENCY).
DO i = il1 , iu1 , 2
Rz(j) = Rz(j) - Amat(i)*Duals(Imat(i)) - Amat(i+1)*Duals(Imat(i+1))
ENDDO
IF ( ihi>Lmx-2 ) THEN
ipage = ipage + 1
key = 1
CALL DPRWPG(key,ipage,lpg,Amat,Imat)
il1 = Nvars + 5
ihi = ihi - lpg
GOTO 1500
ENDIF
ENDIF
pagepl = ihi==(Lmx-2)
!
!     NONBASIC DEPENDENT VARIABLES (COLUMNS DEFINED IMPLICITLY)
Rz(j) = Rz(j)*Csc(j)
!
1600 rcost = Rz(j)
IF ( MOD(Ibb(j),2)==0 ) rcost = -rcost
IF ( Ind(j)==3 ) THEN
IF ( Bu(j)==Bl(j) ) rcost = zero
ENDIF
IF ( Ind(j)==4 ) rcost = -ABS(rcost)
cnorm = one
IF ( j<=Nvars ) cnorm = Colnrm(j)
IF ( rcost+Erdnrm*Dulnrm*cnorm<zero ) nnegrc = nnegrc + 1
j = MOD(j,Mrelas+Nvars) + 1
IF ( nnegrc<Npp.AND.j/=Jstrt ) GOTO 1400
Jstrt = j
1700 RETURN
!     PROCEDURE (TRANSLATE RIGHT HAND SIDE)
!
!     PERFORM THE TRANSLATION ON THE RIGHT-HAND SIDE.
1800 IF ( ibas>Nvars ) THEN
i = ibas - Nvars
IF ( Ind(ibas)/=2 ) THEN
Rhs(i) = Rhs(i) + scalr
ELSE
Rhs(i) = Rhs(i) - scalr
ENDIF
ELSE
i = 0
DO
CALL DPNNZR(i,aij,iplace,Amat,Imat,ibas)
IF ( i<=0 ) EXIT
Rhs(i) = Rhs(i) - scalr*aij*Csc(ibas)
ENDDO
ENDIF
Rhsnrm = MAX(Rhsnrm,DASUM(Mrelas,Rhs,1))
SELECT CASE(npr001)
CASE(100)
GOTO 100
CASE(300)
GOTO 300
CASE(400)
GOTO 400
CASE(800)
GOTO 800
END SELECT
!     PROCEDURE (COMPUTE NEW DUALS)
!
!     SOLVE FOR DUAL VARIABLES. FIRST COPY COSTS INTO DUALS(*).
1900 i = 1
n20121 = Mrelas
DO WHILE ( (n20121-i)>=0 )
j = Ibasis(i)
IF ( j>Nvars ) THEN
Duals(i) = Xlamda*Primal(i+Nvars)
ELSE
Duals(i) = Costsc*Costs(j)*Csc(j) + Xlamda*Primal(i+Nvars)
ENDIF
i = i + 1
ENDDO
!
trans = .TRUE.
CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Duals,trans)
Dulnrm = DASUM(Mrelas,Duals,1)
SELECT CASE(npr003)
CASE(1100)
GOTO 1100
CASE(1300)
GOTO 1300
CASE(1700)
GOTO 1700
END SELECT
END SUBROUTINE DPLPMU
