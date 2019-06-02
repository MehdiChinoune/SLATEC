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
  USE DML, ONLY : c_com, inhomo_com, nfc_com, px_com, pwcnd_com, tnd_com, x_com, &
    xend_com, xot_com, knswot_com, lotjp_com, mnswot_com, nswot_com, tol_com, &
    nps_com, nfcc_com
  !
  INTEGER :: Iflag, Ncomp, Niv, Ip(:)
  REAL(8) :: P(:), S(:), Stowa(:), W(:), Y(:,:), Yhp(:,:), Yp(:)
  INTEGER :: ijk, j, k, kk, l, mflag, nfcp
  REAL(8) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm
  !
  !- *********************************************************************
  !     BEGIN BLOCK PERMITTING ...EXITS TO 210
  !        BEGIN BLOCK PERMITTING ...EXITS TO 10
  !* FIRST EXECUTABLE STATEMENT  DREORT
  nfcp = nfc_com + 1
  !
  !           CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
  !
  !        ...EXIT
  IF ( Iflag==1 ) THEN
    knswot_com = knswot_com + 1
    !        ...EXIT
    IF ( knswot_com<nswot_com ) THEN
      !     ......EXIT
      IF ( (xend_com-x_com)*(x_com-xot_com)<0.0D0 ) RETURN
    END IF
  END IF
  CALL DSTOR1(Y(:,1),Yhp(:,1),Yp,Yhp(:,nfcp),1,0,0)
  !
  !        ***************************************************************
  !
  !        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
  !        AND PARTICULAR SOLUTION YP.
  !
  Niv = nfc_com
  CALL DMGSBV(Ncomp,nfc_com,Y,Ncomp,Niv,mflag,S,P,Ip,inhomo_com,Yp,W,wcnd)
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
      IF ( wcnd>=50.0D0*tol_com ) THEN
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
        nswot_com = knswot_com
        knswot_com = 0
        lotjp_com = 0
        wcnd = LOG10(wcnd)
        IF ( wcnd>tnd_com+3.0D0 ) nswot_com = 2*nswot_com
        IF ( wcnd<pwcnd_com ) THEN
          dx = x_com - px_com
          dnd = pwcnd_com - wcnd
          IF ( dnd>=4 ) nswot_com = nswot_com/2
          dndt = wcnd - tnd_com
          IF ( ABS(dx*dndt)<=dnd*ABS(xend_com-x_com) ) THEN
            xot_com = x_com + dx*dndt/dnd
            nswot_com = MIN(mnswot_com,nswot_com)
            pwcnd_com = wcnd
            !           ......EXIT
            px_com = x_com
          ELSE
            xot_com = xend_com
            nswot_com = MIN(mnswot_com,nswot_com)
            pwcnd_com = wcnd
            px_com = x_com
          END IF
        ELSE
          xot_com = xend_com
          nswot_com = MIN(mnswot_com,nswot_com)
          pwcnd_com = wcnd
          px_com = x_com
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
    50  nswot_com = 1
    knswot_com = 0
    lotjp_com = 1
    kk = 1
    l = 1
    DO k = 1, nfcc_com
      !                 BEGIN BLOCK PERMITTING ...EXITS TO 140
      srp = SQRT(P(kk))
      IF ( inhomo_com==1 ) W(k) = srp*W(k)
      vnorm = 1.0D0/srp
      P(kk) = vnorm
      kk = kk + nfcc_com + 1 - k
      IF ( nfc_com/=nfcc_com ) THEN
        !                 ......EXIT
        IF ( l/=k/2 ) CYCLE
      END IF
      DO j = 1, Ncomp
        Y(j,l) = Y(j,l)*vnorm
      END DO
      l = l + 1
    END DO
    !
    IF ( inhomo_com==1.AND.nps_com/=1 ) THEN
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
      DO j = 1, nfcc_com
        W(j) = c_com*W(j)
      END DO
    END IF
    !
    IF ( Iflag==1 ) CALL DSTWAY(Y(:,1),Yp,Yhp(:,1),0,Stowa)
    Iflag = 0
    !           BEGIN BLOCK PERMITTING ...EXITS TO 40
  ELSEIF ( Iflag==2 ) THEN
    Iflag = 30
  ELSEIF ( nswot_com<=1.AND.lotjp_com/=0 ) THEN
    Iflag = 30
  ELSE
    !
    !                    RETRIEVE DATA FOR A RESTART AT LAST
    !                    ORTHONORMALIZATION POINT
    !
    CALL DSTWAY(Y(:,1),Yp,Yhp(:,1),1,Stowa)
    lotjp_com = 1
    nswot_com = 1
    knswot_com = 0
    mnswot_com = mnswot_com/2
    tnd_com = tnd_com + 1.0D0
    !           .........EXIT
    Iflag = 10
  END IF
  RETURN
END SUBROUTINE DREORT
