!DECK DPLPCE
SUBROUTINE DPLPCE(Mrelas,Nvars,Lmx,Lbm,Itlp,Itbrc,Ibasis,Imat,Ibrc,Ipr,&
    Iwr,Ind,Ibb,Erdnrm,Eps,Tune,Gg,Amat,Basmat,Csc,Wr,Ww,&
    Primal,Erd,Erp,Singlr,Redbas)
  IMPLICIT NONE
  INTEGER i, IDLOC, ihi, il1, ilow, ipage, Itbrc, Itlp, iu1, j, &
    key, l, Lbm, Lmx, lpg, Mrelas, n20002, n20012, n20016, &
    n20023
  INTEGER n20047, n20057, n20061, Nvars
  !***BEGIN PROLOGUE  DPLPCE
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DSPLP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (SPLPCE-S, DPLPCE-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
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
  !***SEE ALSO  DSPLP
  !***ROUTINES CALLED  DASUM, DCOPY, DPRWPG, IDLOC, LA05BD
  !***REVISION HISTORY  (YYMMDD)
  !   811215  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890605  Removed unreferenced labels.  (WRB)
  !   890606  Changed references from IPLOC to IDLOC.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DPLPCE
  INTEGER Ibasis(*), Imat(*), Ibrc(Lbm,2), Ipr(*), Iwr(*), Ind(*), &
    Ibb(*)
  REAL(8) :: Amat(*), Basmat(*), Csc(*), Wr(*), Ww(*), Primal(*)&
    , Erd(*), Erp(*), Eps, Erdnrm, factor, Gg, one, &
    zero, ten, Tune
  REAL(8) :: DASUM
  LOGICAL Singlr, Redbas, trans, pagepl
  !***FIRST EXECUTABLE STATEMENT  DPLPCE
  zero = 0.D0
  one = 1.D0
  ten = 10.D0
  lpg = Lmx - (Nvars+4)
  Singlr = .FALSE.
  factor = 0.01
  !
  !     COPY COLSUMS IN WW(*), AND SOLVE TRANSPOSED SYSTEM.
  i = 1
  n20002 = Mrelas
  DO WHILE ( (n20002-i)>=0 )
    j = Ibasis(i)
    IF ( j<=Nvars ) THEN
      Ww(i) = Primal(j)
    ELSEIF ( Ind(j)/=2 ) THEN
      Ww(i) = -one
    ELSE
      Ww(i) = one
    ENDIF
    i = i + 1
  ENDDO
  !
  !     PERTURB RIGHT-SIDE IN UNITS OF LAST BITS TO BETTER REFLECT
  !     ERRORS IN THE CHECK SUM SOLNS.
  i = 1
  n20012 = Mrelas
  DO WHILE ( (n20012-i)>=0 )
    Ww(i) = Ww(i) + ten*Eps*Ww(i)
    i = i + 1
  ENDDO
  trans = .TRUE.
  CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
  i = 1
  n20016 = Mrelas
  DO WHILE ( (n20016-i)>=0 )
    Erd(i) = MAX(ABS(Ww(i)-one),Eps)*Tune
    !
    !     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS .GT. FACTOR.
    !     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED.
    Singlr = Singlr .OR. (Erd(i)>=factor)
    i = i + 1
  ENDDO
  Erdnrm = DASUM(Mrelas,Erd,1)
  !
  !     RECALCULATE ROW CHECK SUMS EVERY ITBRC ITERATIONS OR WHEN
  !     A REDECOMPOSITION HAS OCCURRED.
  IF ( MOD(Itlp,Itbrc)==0.OR.Redbas ) THEN
    !
    !     COMPUTE ROW SUMS, STORE IN WW(*), SOLVE PRIMAL SYSTEM.
    Ww(1) = zero
    CALL DCOPY(Mrelas,Ww,0,Ww,1)
    pagepl = .TRUE.
    j = 1
    n20023 = Nvars
    DO WHILE ( (n20023-j)>=0 )
      IF ( Ibb(j)<zero ) THEN
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
          DO i = il1, iu1
            Ww(Imat(i)) = Ww(Imat(i)) + Amat(i)*Csc(j)
          ENDDO
          IF ( ihi<=Lmx-2 ) EXIT
          ipage = ipage + 1
          key = 1
          CALL DPRWPG(key,ipage,lpg,Amat,Imat)
          il1 = Nvars + 5
          ihi = ihi - lpg
        ENDDO
        pagepl = ihi==(Lmx-2)
      ELSE
        !
        !     THE VARIABLE IS NON-BASIC.
        pagepl = .TRUE.
      ENDIF
      j = j + 1
    ENDDO
    l = 1
    n20047 = Mrelas
    DO WHILE ( (n20047-l)>=0 )
      j = Ibasis(l)
      IF ( j>Nvars ) THEN
        i = j - Nvars
        IF ( Ind(j)/=2 ) THEN
          Ww(i) = Ww(i) - one
        ELSE
          Ww(i) = Ww(i) + one
        ENDIF
      ENDIF
      l = l + 1
    ENDDO
    !
    !     PERTURB RIGHT-SIDE IN UNITS OF LAST BIT POSITIONS.
    i = 1
    n20057 = Mrelas
    DO WHILE ( (n20057-i)>=0 )
      Ww(i) = Ww(i) + ten*Eps*Ww(i)
      i = i + 1
    ENDDO
    trans = .FALSE.
    CALL LA05BD(Basmat,Ibrc,Lbm,Mrelas,Ipr,Iwr,Wr,Gg,Ww,trans)
    i = 1
    n20061 = Mrelas
    DO WHILE ( (n20061-i)>=0 )
      Erp(i) = MAX(ABS(Ww(i)-one),Eps)*Tune
      !
      !     SYSTEM BECOMES SINGULAR WHEN ACCURACY OF SOLUTION IS .GT. FACTOR.
      !     THIS VALUE (FACTOR) MIGHT NEED TO BE CHANGED.
      Singlr = Singlr .OR. (Erp(i)>=factor)
      i = i + 1
    ENDDO
  ENDIF
  !
END SUBROUTINE DPLPCE
