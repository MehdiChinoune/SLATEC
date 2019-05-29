!** DPINCW
SUBROUTINE DPINCW(Mrelas,Nvars,Lmx,Lbm,Npp,Jstrt,Imat,Ibrc,Ipr,Iwr,&
    Ind,Ibb,Costsc,Gg,Erdnrm,Dulnrm,Amat,Basmat,Csc,Wr,Ww,&
    Rz,Rg,Costs,Colnrm,Duals,Stpedg)
  !>
  !  Subsidiary to DSPLP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SPINCW-S, DPINCW-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
  !     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
  !
  !     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/,
  !     REAL (12 BLANKS)/DOUBLE PRECISION/,/SCOPY/DCOPY/,/SDOT/DDOT/.
  !
  !     THIS SUBPROGRAM IS PART OF THE DSPLP( ) PACKAGE.
  !     IT IMPLEMENTS THE PROCEDURE (INITIALIZE REDUCED COSTS AND
  !     STEEPEST EDGE WEIGHTS).
  !
  !***
  ! **See also:**  DSPLP
  !***
  ! **Routines called:**  DCOPY, DDOT, DPRWPG, IDLOC, LA05BD

  !* REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  INTEGER :: Jstrt, Lbm, Lmx, Mrelas, Npp, Nvars
  REAL(8) :: Costsc, Erdnrm, Dulnrm, Gg
  LOGICAL :: Stpedg
  INTEGER :: Imat(Lmx), Ibrc(Lbm,2), Ipr(2*Mrelas), Iwr(8*Mrelas), &
    Ind(Nvars+Mrelas), Ibb(Nvars+Mrelas)
  REAL(8) :: Amat(Lmx), Basmat(Lbm), Csc(Nvars), Wr(Mrelas), Ww(Mrelas), &
    Rz(Nvars+Mrelas), Rg(Nvars+Mrelas), Costs(Nvars), Colnrm(Nvars), &
    Duals(Nvars+Mrelas)
  INTEGER :: i, ihi, il1, ilow, ipage, iu1, j, key, lpg, nnegrc
  REAL(8) :: one, rzj, scalr, zero, rcost, cnorm
  LOGICAL :: pagepl, trans
  !* FIRST EXECUTABLE STATEMENT  DPINCW
  lpg = Lmx - (Nvars+4)
  zero = 0.D0
  one = 1.D0
  !
  !     FORM REDUCED COSTS, RZ(*), AND STEEPEST EDGE WEIGHTS, RG(*).
  pagepl = .TRUE.
  Rz(1:Nvars+Mrelas) = zero
  Rg(1:Nvars+Mrelas) = one
  nnegrc = 0
  j = Jstrt
  100 CONTINUE
  IF ( Ibb(j)<=0 ) THEN
    pagepl = .TRUE.
    !
    !     THESE ARE NONBASIC INDEPENDENT VARIABLES. THE COLS. ARE IN SPARSE
    !     MATRIX FORMAT.
  ELSEIF ( j>Nvars ) THEN
    pagepl = .TRUE.
    Ww(1:Mrelas) = zero
    scalr = -one
    IF ( Ind(j)==2 ) scalr = one
    i = j - Nvars
    Rz(j) = -scalr*Duals(i)
    Ww(i) = scalr
    IF ( Stpedg ) THEN
      trans = .FALSE.
      CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
      Rg(j) = NORM2(Ww(1:Mrelas))**2 + one
    END IF
  ELSE
    rzj = Costsc*Costs(j)
    Ww(1:Mrelas) = zero
    IF ( j/=1 ) THEN
      ilow = Imat(j+3) + 1
    ELSE
      ilow = Nvars + 5
    END IF
    IF ( .NOT.(pagepl) ) THEN
      il1 = ihi + 1
    ELSE
      il1 = IDLOC(ilow,Amat,Imat)
      IF ( il1>=Lmx-1 ) THEN
        ilow = ilow + 2
        il1 = IDLOC(ilow,Amat,Imat)
      END IF
      ipage = ABS(Imat(Lmx-1))
    END IF
    ihi = Imat(j+4) - (ilow-il1)
    DO
      iu1 = MIN(Lmx-2,ihi)
      IF ( il1>iu1 ) EXIT
      DO i = il1, iu1
        rzj = rzj - Amat(i)*Duals(Imat(i))
        Ww(Imat(i)) = Amat(i)*Csc(j)
      END DO
      IF ( ihi<=Lmx-2 ) EXIT
      ipage = ipage + 1
      key = 1
      CALL DPRWPG(key,ipage,lpg,Amat,Imat)
      il1 = Nvars + 5
      ihi = ihi - lpg
    END DO
    pagepl = ihi==(Lmx-2)
    Rz(j) = rzj*Csc(j)
    IF ( Stpedg ) THEN
      trans = .FALSE.
      CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
      !
      !     THESE ARE NONBASIC DEPENDENT VARIABLES. THE COLS. ARE IMPLICITLY
      !     DEFINED.
      Rg(j) = NORM2(Ww(1:Mrelas))**2 + one
    END IF
  END IF
  !
  rcost = Rz(j)
  IF ( MOD(Ibb(j),2)==0 ) rcost = -rcost
  IF ( Ind(j)==4 ) rcost = -ABS(rcost)
  cnorm = one
  IF ( j<=Nvars ) cnorm = Colnrm(j)
  IF ( rcost+Erdnrm*Dulnrm*cnorm<zero ) nnegrc = nnegrc + 1
  j = MOD(j,Mrelas+Nvars) + 1
  IF ( nnegrc<Npp.AND.j/=Jstrt ) GOTO 100
  Jstrt = j
END SUBROUTINE DPINCW
