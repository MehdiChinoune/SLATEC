!** SPLPDM
SUBROUTINE SPLPDM(Mrelas,Nvars,Lbm,Nredc,Info,Ibasis,Imat,Ibrc,&
    Ipr,Iwr,Ind,Anorm,Eps,Uu,Gg,Amat,Basmat,Csc,Wr,Singlr,Redbas)
  !> Subsidiary to SPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SPLPDM-S, DPLPDM-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THIS SUBPROGRAM IS FROM THE SPLP( ) PACKAGE.  IT PERFORMS THE
  !     TASK OF DEFINING THE ENTRIES OF THE BASIS MATRIX AND
  !     DECOMPOSING IT USING THE LA05 PACKAGE.
  !     IT IS THE MAIN PART OF THE PROCEDURE (DECOMPOSE BASIS MATRIX).
  !
  !***
  ! **See also:**  SPLP
  !***
  ! **Routines called:**  LA05AS, PNNZRS, SASUM, XERMSG
  !***
  ! COMMON BLOCKS    LA05DS

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890605  Corrected references to XERRWV.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls, changed do-it-yourself
  !           DO loops to DO loops.  (RWC)
  USE LA05DS, ONLY : small_com

  INTEGER, INTENT(IN) :: Lbm, Mrelas, Nvars
  INTEGER, INTENT(INOUT) :: Nredc
  INTEGER, INTENT(OUT) :: Info
  REAL(SP), INTENT(IN) :: Eps
  REAL(SP), INTENT(INOUT) :: Uu
  REAL(SP), INTENT(OUT) :: Anorm, Gg
  LOGICAL, INTENT(OUT) :: Singlr, Redbas
  INTEGER, INTENT(IN) :: Ibasis(Nvars+Mrelas), Imat(:), Ind(Nvars+Mrelas)
  INTEGER, INTENT(INOUT) :: Ipr(2*Mrelas)
  INTEGER, INTENT(OUT) :: Ibrc(Lbm,2), Iwr(8*Mrelas)
  REAL(SP), INTENT(IN) :: Amat(:), Csc(Nvars)
  REAL(SP), INTENT(OUT) :: Basmat(Lbm), Wr(Mrelas)
  INTEGER :: i, iplace, j, k, nzbm
  REAL(SP) :: aij
  CHARACTER(16) :: xern3
  !
  !* FIRST EXECUTABLE STATEMENT  SPLPDM
  !
  !     DEFINE BASIS MATRIX BY COLUMNS FOR SPARSE MATRIX EQUATION SOLVER.
  !     THE LA05AS() SUBPROGRAM REQUIRES THE NONZERO ENTRIES OF THE MATRIX
  !     TOGETHER WITH THE ROW AND COLUMN INDICES.
  !
  nzbm = 0
  !
  !     DEFINE DEPENDENT VARIABLE COLUMNS. THESE ARE
  !     COLS. OF THE IDENTITY MATRIX AND IMPLICITLY GENERATED.
  !
  DO k = 1, Mrelas
    j = Ibasis(k)
    IF( j>Nvars ) THEN
      nzbm = nzbm + 1
      IF( Ind(j)==2 ) THEN
        Basmat(nzbm) = 1._SP
      ELSE
        Basmat(nzbm) = -1._SP
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
        CALL PNNZRS(i,aij,iplace,Amat,Imat,j)
        IF( i>0 ) THEN
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
  CALL LA05AS(Basmat,Ibrc,nzbm,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Uu)
  !
  !     CHECK RETURN VALUE OF ERROR FLAG, GG.
  !
  IF( Gg>=0._SP ) RETURN
  IF( Gg==(-7.) ) THEN
    Info = -28
    ERROR STOP 'SPLPDM : IN SPLP, SHORT ON STORAGE FOR LA05AS. USE PRGOPT(*) TO GIVE MORE.'
  ELSEIF( Gg==(-5.) ) THEN
    Singlr = .TRUE.
  ELSE
    WRITE (xern3,'(1PE15.6)') Gg
    Info = -27
    ERROR STOP 'SPLPDM : IN SPLP, LA05AS RETURNED ERROR FLAG' ! = '//xern3
  END IF

END SUBROUTINE SPLPDM