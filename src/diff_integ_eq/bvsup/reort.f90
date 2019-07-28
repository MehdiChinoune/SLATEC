!** REORT
SUBROUTINE REORT(Ncomp,Y,Yp,Yhp,Niv,W,S,P,Ip,Stowa,Iflag)
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (REORT-S, DREORT-D)
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
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  MGSBV, SDOT, STOR1, STWAY
  !***
  ! COMMON BLOCKS    ML15TO, ML18JR, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY : c_com, inhomo_com, nfc_com, px_com, pwcnd_com, tnd_com, x_com, &
    xend_com, xot_com, knswot_com, lotjp_com, mnswot_com, nswot_com, tol_com, &
    nps_com, nfcc_com
  !
  INTEGER, INTENT(IN) :: Ncomp
  INTEGER, INTENT(INOUT) :: Iflag
  INTEGER, INTENT(OUT) :: Niv, Ip(:)
  REAL(SP), INTENT(INOUT) :: Y(:,:), Yhp(:,:), Yp(:)
  REAL(SP), INTENT(OUT) :: P(:), S(:), Stowa(:), W(:)
  !
  INTEGER :: nfcp,ijk, j, k, kk, l, mflag
  REAL(SP) :: dnd, dndt, dx, srp, vnorm, wcnd, ypnm
  !
  !- *********************************************************************
  !* FIRST EXECUTABLE STATEMENT  REORT
  nfcp = nfc_com + 1
  !
  !     CHECK TO SEE IF ORTHONORMALIZATION TEST IS TO BE PERFORMED
  !
  IF( Iflag==1 ) THEN
    knswot_com = knswot_com + 1
    IF( knswot_com<nswot_com ) THEN
      IF( (xend_com-x_com)*(x_com-xot_com)<0. ) RETURN
    END IF
  END IF
  CALL STOR1(Y(:,1),Yhp(:,1),Yp,Yhp(:,nfcp),1,0,0)
  !
  !     ****************************************
  !
  !     ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
  !     AND PARTICULAR SOLUTION YP.
  !
  Niv = nfc_com
  CALL MGSBV(Ncomp,nfc_com,Y,Ncomp,Niv,mflag,S,P,Ip,inhomo_com,Yp,W,wcnd)
  !
  !     ****************************************
  !
  !  CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
  !
  IF( mflag==0 ) THEN
    !
    !     ****************************************
    !
    IF( Iflag==1 ) THEN
      !
      !     TEST FOR ORTHONORMALIZATION
      !
      IF( wcnd>=50._SP*tol_com ) THEN
        DO ijk = 1, nfcp
          IF( S(ijk)>1.0E+20 ) GOTO 50
        END DO
        !
        !     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
        !     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
        !     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
        !     ARE ADDED FOR SAFETY PURPOSES.
        !
        nswot_com = knswot_com
        knswot_com = 0
        lotjp_com = 0
        wcnd = LOG10(wcnd)
        IF( wcnd>tnd_com+3. ) nswot_com = 2*nswot_com
        IF( wcnd>=pwcnd_com ) THEN
          xot_com = xend_com
        ELSE
          dx = x_com - px_com
          dnd = pwcnd_com - wcnd
          IF( dnd>=4 ) nswot_com = nswot_com/2
          dndt = wcnd - tnd_com
          IF( ABS(dx*dndt)>dnd*ABS(xend_com-x_com) ) THEN
            xot_com = xend_com
          ELSE
            xot_com = x_com + dx*dndt/dnd
          END IF
        END IF
        nswot_com = MIN(mnswot_com,nswot_com)
        pwcnd_com = wcnd
        px_com = x_com
        RETURN
      END IF
    END IF
    !
    !     ****************************************
    !
    !     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
    !     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
    !
    50  nswot_com = 1
    knswot_com = 0
    lotjp_com = 1
    kk = 1
    l = 1
    DO k = 1, nfcc_com
      srp = SQRT(P(kk))
      IF( inhomo_com==1 ) W(k) = srp*W(k)
      vnorm = 1._SP/srp
      P(kk) = vnorm
      kk = kk + nfcc_com + 1 - k
      IF( nfc_com/=nfcc_com ) THEN
        IF( l/=k/2 ) CYCLE
      END IF
      DO j = 1, Ncomp
        Y(j,l) = Y(j,l)*vnorm
      END DO
      l = l + 1
    END DO
    !
    IF( inhomo_com==1 .AND. nps_com/=1 ) THEN
      !
      !     NORMALIZE THE PARTICULAR SOLUTION
      !
      ypnm = NORM2(Yp(1:Ncomp))**2
      IF( ypnm==0._SP ) ypnm = 1._SP
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
    IF( Iflag==1 ) CALL STWAY(Y(:,1),Yp,Yhp(:,1),0,Stowa)
    Iflag = 0
  ELSE
    IF( Iflag/=2 ) THEN
      IF( nswot_com>1 .OR. lotjp_com==0 ) THEN
        !
        !     RETRIEVE DATA FOR A RESTART AT LAST ORTHONORMALIZATION POINT
        !
        CALL STWAY(Y(:,1),Yp,Yhp(:,1),1,Stowa)
        lotjp_com = 1
        nswot_com = 1
        knswot_com = 0
        mnswot_com = mnswot_com/2
        tnd_com = tnd_com + 1._SP
        Iflag = 10
        RETURN
      END IF
    END IF
    Iflag = 30
    RETURN
  END IF
  !
END SUBROUTINE REORT