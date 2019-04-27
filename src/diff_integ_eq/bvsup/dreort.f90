!** DREORT
SUBROUTINE DREORT(Ncomp,Y,Yp,Yhp,Niv,W,S,P,Ip,Stowa,Iflag)
  !>
  !  Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (REORT-S, DREORT-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !   INPUT
  !- ********
  !     Y, YP and YHP = homogeneous solution matrix and particular
  !                     solution vector to be orthonormalized.
  !     IFLAG = 1 --  store YHP into Y and YP, test for
  !                   reorthonormalization, orthonormalize if needed,
  !                   save restart data.
  !             2 --  store YHP into Y and YP, reorthonormalization,
  !                   no restarts.
  !                   (preset orthonormalization mode)
  !             3 --  store YHP into Y and YP, reorthonormalization
  !                   (when INHOMO=3 and X=XEND).
  !- *********************************************************************
  !   OUTPUT
  !- ********
  !     Y, YP = orthonormalized solutions.
  !     NIV = number of independent vectors returned from DMGSBV.
  !     IFLAG = 0 --  reorthonormalization was performed.
  !            10 --  solution process must be restarted at the last
  !                   orthonormalization point.
  !            30 --  solutions are linearly dependent, problem must
  !                   be restarted from the beginning.
  !     W, P, IP = orthonormalization information.
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DDOT, DMGSBV, DSTOR1, DSTWAY
  !***
  ! COMMON BLOCKS    DML15T, DML18J, DML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : C, INHomo, NFC, PX, PWCnd, TND, X, XENd, XOT, KNSwot, LOTjp, &
    MNSwot, NSWot, TOL, NPS, NFCc
  !
  INTEGER Iflag, ijk, Ip(*), j, k, kk, l, mflag, Ncomp, nfcp, Niv
  REAL(8) :: dnd, dndt, dx, P(*), S(*), srp, Stowa(*), vnorm, W(*), wcnd, &
    Y(Ncomp,*), Yhp(Ncomp,*), Yp(*), ypnm
  !
  !- *********************************************************************
  !     BEGIN BLOCK PERMITTING ...EXITS TO 210
  !        BEGIN BLOCK PERMITTING ...EXITS TO 10
  !* FIRST EXECUTABLE STATEMENT  DREORT
  nfcp = NFC + 1
  !
  !           CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
  !
  !        ...EXIT
  IF ( Iflag==1 ) THEN
    KNSwot = KNSwot + 1
    !        ...EXIT
    IF ( KNSwot<NSWot ) THEN
      !     ......EXIT
      IF ( (XENd-X)*(X-XOT)<0.0D0 ) RETURN
    END IF
  END IF
  CALL DSTOR1(Y,Yhp,Yp,Yhp(1,nfcp),1,0,0)
  !
  !        ***************************************************************
  !
  !        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
  !        AND PARTICULAR SOLUTION YP.
  !
  Niv = NFC
  CALL DMGSBV(Ncomp,NFC,Y,Ncomp,Niv,mflag,S,P,Ip,INHomo,Yp,W,wcnd)
  !
  !           ************************************************************
  !
  !        CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
  !
  IF ( mflag==0 ) THEN
    !           BEGIN BLOCK PERMITTING ...EXITS TO 190
    !              BEGIN BLOCK PERMITTING ...EXITS TO 110
    !
    !                 ******************************************************
    !
    !              ...EXIT
    IF ( Iflag==1 ) THEN
      !
      !                 TEST FOR ORTHONORMALIZATION
      !
      !              ...EXIT
      IF ( wcnd>=50.0D0*TOL ) THEN
        DO ijk = 1, nfcp
          !              ......EXIT
          IF ( S(ijk)>1.0D20 ) GOTO 50
        END DO
        !
        !                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
        !                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
        !                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
        !                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
        !
        NSWot = KNSwot
        KNSwot = 0
        LOTjp = 0
        wcnd = LOG10(wcnd)
        IF ( wcnd>TND+3.0D0 ) NSWot = 2*NSWot
        IF ( wcnd<PWCnd ) THEN
          dx = X - PX
          dnd = PWCnd - wcnd
          IF ( dnd>=4 ) NSWot = NSWot/2
          dndt = wcnd - TND
          IF ( ABS(dx*dndt)<=dnd*ABS(XENd-X) ) THEN
            XOT = X + dx*dndt/dnd
            NSWot = MIN(MNSwot,NSWot)
            PWCnd = wcnd
            !           ......EXIT
            PX = X
          ELSE
            XOT = XENd
            NSWot = MIN(MNSwot,NSWot)
            PWCnd = wcnd
            PX = X
          END IF
        ELSE
          XOT = XENd
          NSWot = MIN(MNSwot,NSWot)
          PWCnd = wcnd
          PX = X
        END IF
        RETURN
      END IF
    END IF
    !
    !              *********************************************************
    !
    !              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
    !              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50  NSWot = 1
    KNSwot = 0
    LOTjp = 1
    kk = 1
    l = 1
    DO k = 1, NFCc
      !                 BEGIN BLOCK PERMITTING ...EXITS TO 140
      srp = SQRT(P(kk))
      IF ( INHomo==1 ) W(k) = srp*W(k)
      vnorm = 1.0D0/srp
      P(kk) = vnorm
      kk = kk + NFCc + 1 - k
      IF ( NFC/=NFCc ) THEN
        !                 ......EXIT
        IF ( l/=k/2 ) CYCLE
      END IF
      DO j = 1, Ncomp
        Y(j,l) = Y(j,l)*vnorm
      END DO
      l = l + 1
    END DO
    !
    IF ( INHomo==1.AND.NPS/=1 ) THEN
      !
      !                 NORMALIZE THE PARTICULAR SOLUTION
      !
      ypnm = NORM2(Yp(1:Ncomp))**2
      IF ( ypnm==0.0D0 ) ypnm = 1.0D0
      ypnm = SQRT(ypnm)
      S(nfcp) = ypnm
      DO j = 1, Ncomp
        Yp(j) = Yp(j)/ypnm
      END DO
      DO j = 1, NFCc
        W(j) = C*W(j)
      END DO
    END IF
    !
    IF ( Iflag==1 ) CALL DSTWAY(Y,Yp,Yhp,0,Stowa)
    Iflag = 0
    !           BEGIN BLOCK PERMITTING ...EXITS TO 40
  ELSEIF ( Iflag==2 ) THEN
    Iflag = 30
  ELSEIF ( NSWot<=1.AND.LOTjp/=0 ) THEN
    Iflag = 30
  ELSE
    !
    !                    RETRIEVE DATA FOR A RESTART AT LAST
    !                    ORTHONORMALIZATION POINT
    !
    CALL DSTWAY(Y,Yp,Yhp,1,Stowa)
    LOTjp = 1
    NSWot = 1
    KNSwot = 0
    MNSwot = MNSwot/2
    TND = TND + 1.0D0
    !           .........EXIT
    Iflag = 10
  END IF
  RETURN
END SUBROUTINE DREORT
