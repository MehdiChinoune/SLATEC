!** RKFAB
SUBROUTINE RKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,&
    W,S,Stowa,Work,Iwork,Nfcc)
  !> Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (RKFAB-S, DRKFAB-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !
  !     Subroutine RKFAB integrates the initial value equations using
  !     the variable-step RUNGE-KUTTA-FEHLBERG integration scheme or
  !     the variable-order ADAMS method and orthonormalization
  !     determined by a linear dependence test.
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  BVDER, DEABM, DERKF, REORT, STOR1
  !***
  ! COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE ML, ONLY : c_com, inhomo_com, x_com, xbeg_com, xend_com, xop_com, info_com, &
    kop_com, ae_com, re_com, nopg_com, ndisk_com, ntape_com, neq_com, integ_com, &
    nps_com, numort_com, kkkint_com, lllint_com
  INTEGER :: Ncomp, Nfc, Nfcc, Niv, Ntp, Iflag, Mxnon, Nxpts
  INTEGER :: Iwork(*), Ip(Nfcc,Mxnon+1)
  REAL(SP) :: P(Ntp,Mxnon+1), S(Nfc+1), Stowa(:), U(Ncomp,Nfc,Nxpts), &
    V(Ncomp,Nxpts), W(Nfcc,Mxnon+1), Work(*), Xpts(Nxpts), Yhp(Ncomp,Nfc+1), Z(Mxnon+1)
  INTEGER :: nfcp1, non, idid, ipar(1), j, jflag, jon, kod, kopp
  REAL(SP) :: xxop, ret(1), aet(1)
  !
  !- *********************************************************************
  !  INITIALIZATION OF COUNTERS AND VARIABLES.
  !
  !* FIRST EXECUTABLE STATEMENT  RKFAB
  kod = 1
  non = 1
  x_com = xbeg_com
  jon = 1
  info_com(1) = 0
  info_com(2) = 0
  info_com(3) = 1
  info_com(4) = 1
  Work(1) = xend_com
  ipar = 0
  IF( nopg_com/=0 ) THEN
    info_com(3) = 0
    IF( x_com==Z(1) ) jon = 2
  END IF
  nfcp1 = Nfc + 1
  !
  !- *********************************************************************
  !- ****BEGINNING OF INTEGRATION LOOP AT OUTPUT POINTS.******************
  !- *********************************************************************
  !
  DO kopp = 2, Nxpts
    kop_com = kopp
    !
    50  xop_com = Xpts(kop_com)
    IF( ndisk_com==0 ) kod = kop_com
    !
    !     STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
    !
    100  xxop = xop_com
    IF( nopg_com/=0 ) THEN
      IF( xend_com>xbeg_com .AND. xop_com>Z(jon) ) xxop = Z(jon)
      IF( xend_com<xbeg_com .AND. xop_com<Z(jon) ) xxop = Z(jon)
    END IF
    !
    !- *********************************************************************
    150 CONTINUE
    IF( integ_com==2 ) THEN
      !     DEABM INTEGRATOR
      !
      ret(1) = re_com
      aet(1) = ae_com
      CALL DEABM(BVDER_2,neq_com,x_com,Yhp,xxop,info_com,ret,aet,idid,Work,kkkint_com,Iwork,&
        lllint_com)
    ELSE
      !     DERKF INTEGRATOR
      !
      ret(1) = re_com
      aet(1) = ae_com
      CALL DERKF(BVDER_2,neq_com,x_com,Yhp,xxop,info_com,ret,aet,idid,Work,kkkint_com,Iwork,&
        lllint_com)
    END IF
    IF( idid>=1 ) THEN
      !
      !- *********************************************************************
      !     GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR ORTHONORMALIZATION
      !     (TEMPORARILY USING U AND V IN THE TEST)
      !
      IF( nopg_com==0 ) THEN
        jflag = 1
        IF( inhomo_com==3 .AND. x_com==xend_com ) jflag = 3
      ELSE
        IF( xxop/=Z(jon) ) GOTO 200
        jflag = 2
      END IF
      !
      IF( ndisk_com==0 ) non = numort_com + 1
      CALL REORT(Ncomp,U(:,:,kod),V(:,kod),Yhp,Niv,W(:,non),S,P(:,non),&
        Ip(:,non),Stowa,jflag)
      !
      IF( jflag/=30 ) THEN
        !
        IF( jflag==10 ) GOTO 50
        !
        IF( jflag==0 ) THEN
          !
          !- *********************************************************************
          !     STORE ORTHONORMALIZED VECTORS INTO SOLUTION VECTORS.
          !
          IF( numort_com>=Mxnon ) THEN
            IF( x_com/=xend_com ) THEN
              Iflag = 13
              RETURN
            END IF
          END IF
          !
          numort_com = numort_com + 1
          CALL STOR1(Yhp(:,1),U(:,1,kod),Yhp(:,nfcp1),V(:,kod),1,ndisk_com,ntape_com)
          !
          !- *********************************************************************
          !     STORE ORTHONORMALIZATION INFORMATION, INITIALIZE
          !     INTEGRATION FLAG, AND CONTINUE INTEGRATION TO THE NEXT
          !     ORTHONORMALIZATION POINT OR OUTPUT POINT.
          !
          Z(numort_com) = x_com
          IF( inhomo_com==1 .AND. nps_com==0 ) c_com = S(nfcp1)*c_com
          IF( ndisk_com/=0 ) THEN
            IF( inhomo_com==1 ) WRITE (ntape_com) (W(j,1),j=1,Nfcc)
            WRITE (ntape_com) (Ip(j,1),j=1,Nfcc), (P(j,1),j=1,Ntp)
          END IF
          info_com(1) = 0
          jon = jon + 1
          IF( nopg_com==1 .AND. x_com/=xop_com ) GOTO 100
        END IF
      ELSE
        Iflag = 30
        RETURN
      END IF
    ELSE
      info_com(1) = 1
      IF( idid==-1 ) GOTO 150
      Iflag = 20 - idid
      RETURN
    END IF
    !
    !- *********************************************************************
    !     CONTINUE INTEGRATION IF WE ARE NOT AT AN OUTPUT POINT.
    !
    200  IF( idid==1 ) GOTO 150
    !
    !     STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
    !     SOLUTION IN V AT THE OUTPUT POINTS.
    !
    CALL STOR1(U(:,1,kod),Yhp(:,1),V(:,kod),Yhp(:,nfcp1),0,ndisk_com,ntape_com)
  END DO
  !- *********************************************************************
  !- *********************************************************************
  !
  Iflag = 0

CONTAINS
  SUBROUTINE BVDER_2(X,Y,Yp)
    REAL(SP) :: X
    REAL(SP) :: Y(:), Yp(:)
    REAL(SP) :: g(SIZE(Y))

    g = 0._SP
    CALL BVDER(X,Y,Yp,G)
  END SUBROUTINE BVDER_2
END SUBROUTINE RKFAB
