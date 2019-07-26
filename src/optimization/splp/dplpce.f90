!** DPLPCE
PURE SUBROUTINE DPLPCE(Mrelas,Nvars,Lmx,Lbm,Itlp,Itbrc,Ibasis,Imat,Ibrc,Ipr,&
    Iwr,Ind,Ibb,Erdnrm,Eps,Tune,Gg,Amat,Basmat,Csc,Wr,Ww,Primal,Erd,Erp,Singlr,Redbas)
  !> Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SPLPCE-S, DPLPCE-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
  !     /REAL (12 BLANKS)/DOUBLE PRECISION/,
  !     /SASUM/DASUM/,/DCOPY/,DCOPY/.
  !
  !     REVISED 811219-1630
  !     REVISED YYMMDD-HHMM
  !
  !     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT CALCULATES
  !     THE APPROXIMATE ERROR IN THE PRIMAL AND DUAL SYSTEMS.  IT IS
  !     THE MAIN PART OF THE PROCEDURE (COMPUTE ERROR IN DUAL AND PRIMAL
  !     SYSTEMS).
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DASUM, DCOPY, DPRWPG, IDLOC, LA05BD

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  INTEGER, INTENT(IN) :: Itbrc, Itlp, Lbm, Lmx, Mrelas, Nvars
  REAL(DP), INTENT(IN) :: Eps, Gg, Tune
  REAL(DP), INTENT(OUT) :: Erdnrm
  LOGICAL, INTENT(IN) :: Redbas
  LOGICAL, INTENT(OUT) :: Singlr
  INTEGER, INTENT(IN) :: Ibasis(Nvars+Mrelas), Ibrc(Lbm,2), Imat(Lmx), &
    Iwr(8*Mrelas), Ind(Nvars+Mrelas), Ibb(Nvars+Mrelas)
  INTEGER, INTENT(INOUT) :: Ipr(2*Mrelas)
  REAL(DP), INTENT(IN) :: Amat(Lmx), Basmat(Lbm), Csc(Nvars), Primal(Nvars+Mrelas)
  REAL(DP), INTENT(OUT) :: Erd(Mrelas), Erp(Mrelas), Wr(Mrelas), Ww(Mrelas)
  INTEGER :: i, ihi, il1, ilow, ipage, iu1, j, l, lpg, n20002, n20012, &
    n20016, n20023, n20047, n20057, n20061
  REAL(DP) :: factor
  LOGICAL :: trans, pagepl
  !* FIRST EXECUTABLE STATEMENT  DPLPCE
  lpg = Lmx - (Nvars+4)
  Singlr = .FALSE.
  factor = 0.01_DP
  !
  !     COPY COLSUMS IN WW(*), AND SOLVE TRANSPOSED SYSTEM.
  i = 1
  n20002 = Mrelas
  DO WHILE( (n20002-i)>=0 )
    j = Ibasis(i)
    IF( j<=Nvars ) THEN
      Ww(i) = Primal(j)
    ELSEIF( Ind(j)/=2 ) THEN
      Ww(i) = -1._DP
    ELSE
      Ww(i) = 1._DP
    END IF
    i = i + 1
  END DO
  !
  !     PERTURB RIGHT-SIDE IN UNITS OF LAST BITS TO BETTER REFLECT
  !     ERRORS IN THE CHECK SUM SOLNS.
  i = 1
  n20012 = Mrelas
  DO WHILE( (n20012-i)>=0 )
    Ww(i) = Ww(i) + 10._DP*Eps*Ww(i)
    i = i + 1
  END DO
  trans = .TRUE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
  i = 1
  n20016 = Mrelas
  DO WHILE( (n20016-i)>=0 )
    Erd(i) = MAX(ABS(Ww(i)-1._DP),Eps)*Tune
    !
    !     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS > FACTOR.
    !     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED.
    Singlr = Singlr .OR. (Erd(i)>=factor)
    i = i + 1
  END DO
  Erdnrm = SUM(ABS(Erd(1:Mrelas)))
  !
  !     RECALCULATE ROW CHECK SUMS EVERY ITBRC ITERATIONS OR WHEN
  !     A REDECOMPOSITION HAS OCCURRED.
  IF( MOD(Itlp,Itbrc)==0 .OR. Redbas ) THEN
    !
    !     COMPUTE ROW SUMS, STORE IN WW(*), SOLVE PRIMAL SYSTEM.
    Ww(1:Mrelas) = 0._DP
    pagepl = .TRUE.
    j = 1
    n20023 = Nvars
    DO WHILE( (n20023-j)>=0 )
      IF( Ibb(j)<0._DP ) THEN
        IF( j/=1 ) THEN
          ilow = Imat(j+3) + 1
        ELSE
          ilow = Nvars + 5
        END IF
        IF( .NOT. (pagepl) ) THEN
          il1 = ihi + 1
        ELSE
          il1 = IDLOC(ilow,Imat)
          IF( il1>=Lmx-1 ) THEN
            ilow = ilow + 2
            il1 = IDLOC(ilow,Imat)
          END IF
          ipage = ABS(Imat(Lmx-1))
        END IF
        ihi = Imat(j+4) - (ilow-il1)
        DO
          iu1 = MIN(Lmx-2,ihi)
          IF( il1>iu1 ) EXIT
          DO i = il1, iu1
            Ww(Imat(i)) = Ww(Imat(i)) + Amat(i)*Csc(j)
          END DO
          IF( ihi<=Lmx-2 ) EXIT
          il1 = Nvars + 5
          ihi = ihi - lpg
        END DO
        pagepl = ihi==(Lmx-2)
      ELSE
        !
        !     THE VARIABLE IS NON-BASIC.
        pagepl = .TRUE.
      END IF
      j = j + 1
    END DO
    l = 1
    n20047 = Mrelas
    DO WHILE( (n20047-l)>=0 )
      j = Ibasis(l)
      IF( j>Nvars ) THEN
        i = j - Nvars
        IF( Ind(j)/=2 ) THEN
          Ww(i) = Ww(i) - 1._DP
        ELSE
          Ww(i) = Ww(i) + 1._DP
        END IF
      END IF
      l = l + 1
    END DO
    !
    !     PERTURB RIGHT-SIDE IN UNITS OF LAST BIT POSITIONS.
    i = 1
    n20057 = Mrelas
    DO WHILE( (n20057-i)>=0 )
      Ww(i) = Ww(i) + 10._DP*Eps*Ww(i)
      i = i + 1
    END DO
    trans = .FALSE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
    i = 1
    n20061 = Mrelas
    DO WHILE( (n20061-i)>=0 )
      Erp(i) = MAX(ABS(Ww(i)-1._DP),Eps)*Tune
      !
      !     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS > FACTOR.
      !     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED.
      Singlr = Singlr .OR. (Erp(i)>=factor)
      i = i + 1
    END DO
  END IF
  !
END SUBROUTINE DPLPCE