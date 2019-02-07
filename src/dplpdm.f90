!*==DPLPDM.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DPLPDM
SUBROUTINE DPLPDM(Mrelas,Nvars,Lmx,Lbm,Nredc,Info,Iopt,Ibasis,Imat,Ibrc,&
    Ipr,Iwr,Ind,Ibb,Anorm,Eps,Uu,Gg,Amat,Basmat,Csc,Wr,&
    Singlr,Redbas)
  IMPLICIT NONE
  !*--DPLPDM7
  !*** Start of declarations inserted by SPAG
  INTEGER i , Info , Iopt , iplace , j , k , Lbm , LCOl , LENl , LENu ,&
    Lmx , LP , LROw , Mrelas , NCP , Nredc , Nvars , nzbm
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DPLPDM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SPLPDM-S, DPLPDM-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE
  !     TASK OF DEFINING THE ENTRIES OF THE BASIS MATRIX AND
  !     DECOMPOSING IT USING THE LA05 PACKAGE.
  !     IT IS THE MAIN PART OF THE PROCEDURE (DECOMPOSE BASIS MATRIX).
  !
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  DASUM, DPNNZR, LA05AD, XERMSG
  !***COMMON BLOCKS    LA05DD
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Added DASUM to list of DOUBLE PRECISION variables.
  !   890605  Removed unreferenced labels.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls, convert do-it-yourself
  !           DO loops to DO loops.  (RWC)
  !***END PROLOGUE  DPLPDM
  INTEGER Ibasis(*) , Imat(*) , Ibrc(Lbm,2) , Ipr(*) , Iwr(*) , Ind(*) ,&
    Ibb(*)
  REAL(8) :: aij , Amat(*) , Basmat(*) , Csc(*) , Wr(*) , Anorm ,&
    DASUM , Eps , Gg , one , SMAll , Uu , zero
  LOGICAL Singlr , Redbas
  CHARACTER(16) :: xern3
  !
  !     COMMON BLOCK USED BY LA05 () PACKAGE..
  COMMON /LA05DD/ SMAll , LP , LENl , LENu , NCP , LROw , LCOl
  !
  !***FIRST EXECUTABLE STATEMENT  DPLPDM
  zero = 0.D0
  one = 1.D0
  !
  !     DEFINE BASIS MATRIX BY COLUMNS FOR SPARSE MATRIX EQUATION SOLVER.
  !     THE LA05AD() SUBPROGRAM REQUIRES THE NONZERO ENTRIES OF THE MATRIX
  !     TOGETHER WITH THE ROW AND COLUMN INDICES.
  !
  nzbm = 0
  !
  !     DEFINE DEPENDENT VARIABLE COLUMNS. THESE ARE
  !     COLS. OF THE IDENTITY MATRIX AND IMPLICITLY GENERATED.
  !
  DO k = 1 , Mrelas
    j = Ibasis(k)
    IF ( j>Nvars ) THEN
      nzbm = nzbm + 1
      IF ( Ind(j)==2 ) THEN
        Basmat(nzbm) = one
      ELSE
        Basmat(nzbm) = -one
      ENDIF
      Ibrc(nzbm,1) = j - Nvars
      Ibrc(nzbm,2) = k
    ELSE
      !
      !           DEFINE THE INDEP. VARIABLE COLS.  THIS REQUIRES RETRIEVING
      !           THE COLS. FROM THE SPARSE MATRIX DATA STRUCTURE.
      !
      i = 0
      DO
        CALL DPNNZR(i,aij,iplace,Amat,Imat,j)
        IF ( i>0 ) THEN
          nzbm = nzbm + 1
          Basmat(nzbm) = aij*Csc(j)
          Ibrc(nzbm,1) = i
          Ibrc(nzbm,2) = k
          CYCLE
        ENDIF
        EXIT
      ENDDO
    ENDIF
  ENDDO
  !
  Singlr = .FALSE.
  !
  !     RECOMPUTE MATRIX NORM USING CRUDE NORM  =  SUM OF MAGNITUDES.
  !
  Anorm = DASUM(nzbm,Basmat,1)
  SMAll = Eps*Anorm
  !
  !     GET AN L-U FACTORIZATION OF THE BASIS MATRIX.
  !
  Nredc = Nredc + 1
  Redbas = .TRUE.
  CALL LA05AD(Basmat,Ibrc,nzbm,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Uu)
  !
  !     CHECK RETURN VALUE OF ERROR FLAG, GG.
  !
  IF ( Gg>=zero ) RETURN
  IF ( Gg==(-7.) ) THEN
    CALL XERMSG('SLATEC','DPLPDM','IN DSPLP, SHORT ON STORAGE FOR LA05AD.  '&
      //'USE PRGOPT(*) TO GIVE MORE.',28,Iopt)
    Info = -28
  ELSEIF ( Gg==(-5.) ) THEN
    Singlr = .TRUE.
  ELSE
    WRITE (xern3,'(1PE15.6)') Gg
    CALL XERMSG('SLATEC','DPLPDM','IN DSPLP, LA05AD RETURNED ERROR FLAG = '&
      //xern3,27,Iopt)
    Info = -27
  ENDIF
END SUBROUTINE DPLPDM
