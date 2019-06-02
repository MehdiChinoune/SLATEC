!** DPLPDM
SUBROUTINE DPLPDM(Mrelas,Nvars,Lbm,Nredc,Info,Iopt,Ibasis,Imat,Ibrc,&
    Ipr,Iwr,Ind,Anorm,Eps,Uu,Gg,Amat,Basmat,Csc,Wr,Singlr,Redbas)
  !>
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SPLPDM-S, DPLPDM-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE
  !     TASK OF DEFINING THE ENTRIES OF THE BASIS MATRIX AND
  !     DECOMPOSING IT USING THE LA05 PACKAGE.
  !     IT IS THE MAIN PART OF THE PROCEDURE (DECOMPOSE BASIS MATRIX).
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DASUM, DPNNZR, LA05AD, XERMSG
  !***
  ! COMMON BLOCKS    LA05DD

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Added DASUM to list of DOUBLE PRECISION variables.
  !   890605  Removed unreferenced labels.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls, convert do-it-yourself
  !           DO loops to DO loops.  (RWC)
  USE LA05DD, ONLY : small_com
  USE service, ONLY : XERMSG
  INTEGER :: Info, Iopt, Lbm, Mrelas, Nredc, Nvars
  REAL(8) :: Anorm, Eps, Gg, Uu
  LOGICAL :: Singlr, Redbas
  INTEGER :: Ibasis(Nvars+Mrelas), Imat(:), Ibrc(Lbm,2), Ipr(2*Mrelas), &
    Iwr(8*Mrelas), Ind(Nvars+Mrelas)
  REAL(8) :: Amat(:), Basmat(Lbm), Csc(Nvars), Wr(Mrelas)
  INTEGER :: i, iplace, j, k, nzbm
  REAL(8) :: aij, one, zero
  CHARACTER(16) :: xern3
  !
  !* FIRST EXECUTABLE STATEMENT  DPLPDM
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
  DO k = 1, Mrelas
    j = Ibasis(k)
    IF ( j>Nvars ) THEN
      nzbm = nzbm + 1
      IF ( Ind(j)==2 ) THEN
        Basmat(nzbm) = one
      ELSE
        Basmat(nzbm) = -one
      END IF
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
        END IF
        EXIT
      END DO
    END IF
  END DO
  !
  Singlr = .FALSE.
  !
  !     RECOMPUTE MATRIX NORM USING CRUDE NORM  =  SUM OF MAGNITUDES.
  !
  Anorm = SUM(ABS(Basmat(1:nzbm)))
  small_com = Eps*Anorm
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
    CALL XERMSG('SLATEC','DPLPDM','IN DSPLP, SHORT ON STORAGE FOR LA05AD.  USE PRGOPT(*) TO GIVE MORE.',28,Iopt)
    Info = -28
  ELSEIF ( Gg==(-5.) ) THEN
    Singlr = .TRUE.
  ELSE
    WRITE (xern3,'(1PE15.6)') Gg
    CALL XERMSG('SLATEC','DPLPDM','IN DSPLP, LA05AD RETURNED ERROR FLAG = '&
      //xern3,27,Iopt)
    Info = -27
  END IF
END SUBROUTINE DPLPDM
