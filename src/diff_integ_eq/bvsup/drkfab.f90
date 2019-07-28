!** DRKFAB
SUBROUTINE DRKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,&
    W,S,Stowa,Work,Iwork,Nfcc)
  !> Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (RKFAB-S, DRKFAB-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !- *********************************************************************
  !
  !     Subroutine DRKFAB integrates the initial value equations using
  !     the variable-step Runge-Kutta-Fehlberg integration scheme or
  !     the variable-order Adams method and orthonormalization
  !     determined by a linear dependence test.
  !
  !- *********************************************************************
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DBVDER, DDEABM, DDERKF, DREORT, DSTOR1
  !***
  ! COMMON BLOCKS    DML15T, DML17B, DML18J, DML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : c_com, inhomo_com, kkkint_com, lllint_com, x_com, xbeg_com, &
    xend_com, xop_com, info_com, kop_com, ae_com, re_com, nopg_com, ndisk_com, &
    ntape_com, neq_com, integ_com, nps_com, numort_com
  !
  INTEGER, INTENT(IN) :: Mxnon, Ncomp, Nfc, Nfcc, Ntp, Nxpts
  INTEGER, INTENT(INOUT) :: Iflag, Niv
  INTEGER, INTENT(INOUT) :: Ip(Nfcc,Mxnon+1), Iwork(*)
  REAL(DP), INTENT(IN) :: Xpts(Nxpts)
  REAL(DP), INTENT(INOUT) :: P(Ntp,Mxnon+1), S(Nfc+1), Stowa(:), U(Ncomp,Nfc,Nxpts), &
    V(Ncomp,Nxpts), W(Nfcc,Mxnon+1), Work(*), Yhp(Ncomp,Nfc+1), Z(Mxnon+1)
  !
  INTEGER :: idid, ipar(1), j, jflag, jon, kod, kopp, nfcp1, non
  REAL(DP) :: xxop, ret(1), aet(1)
  !
  !      *****************************************************************
  !       INITIALIZATION OF COUNTERS AND VARIABLES.
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 220
  !        BEGIN BLOCK PERMITTING ...EXITS TO 10
  !* FIRST EXECUTABLE STATEMENT  DRKFAB
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
  !        ...EXIT
  IF( nopg_com/=0 ) THEN
    info_com(3) = 0
    IF( x_com==Z(1) ) jon = 2
  END IF
  nfcp1 = Nfc + 1
  !
  !        ***************************************************************
  !        *****BEGINNING OF INTEGRATION LOOP AT OUTPUT
  !        POINTS.******************
  !        ***************************************************************
  !
  DO kopp = 2, Nxpts
    kop_com = kopp
    xop_com = Xpts(kop_com)
    IF( ndisk_com==0 ) kod = kop_com
    !
    !
    !              STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
    !
    !              BEGIN BLOCK PERMITTING ...EXITS TO 190
    !                 BEGIN BLOCK PERMITTING ...EXITS TO 30
    50  xxop = xop_com
    !                 ...EXIT
    IF( nopg_com/=0 ) THEN
      IF( xend_com>xbeg_com .AND. xop_com>Z(jon) ) xxop = Z(jon)
      IF( xend_com<xbeg_com .AND. xop_com<Z(jon) ) xxop = Z(jon)
    END IF
    DO
      !
      !                 ******************************************************
      !                    BEGIN BLOCK PERMITTING ...EXITS TO 170
      IF( integ_com==2 ) THEN
        !                       DDEABM INTEGRATOR
        !
        ret(1) = re_com
        aet(1) = ae_com
        CALL DDEABM(DBVDER_2,neq_com,x_com,Yhp,xxop,info_com,ret,aet,idid,Work,&
          kkkint_com,Iwork,lllint_com)
      ELSE
        !                       DDERKF INTEGRATOR
        !
        ret(1) = re_com
        aet(1) = ae_com
        CALL DDERKF(DBVDER_2,neq_com,x_com,Yhp,xxop,info_com,ret,aet,idid,Work,&
          kkkint_com,Iwork,lllint_com)
      END IF
      IF( idid>=1 ) THEN
        !
        !                       ************************************************
        !                           GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR
        !                           ORTHONORMALIZATION (TEMPORARILY USING U AND
        !                           V IN THE TEST)
        !
        IF( nopg_com==0 ) THEN
          jflag = 1
          IF( inhomo_com==3 .AND. x_com==xend_com ) jflag = 3
        ELSEIF( xxop==Z(jon) ) THEN
          jflag = 2
        ELSE
          !
          !                             ******************************************
          !                                 CONTINUE INTEGRATION IF WE ARE NOT AT
          !                                 AN OUTPUT POINT.
          !
          !           ..................EXIT
          !                    .........EXIT
          IF( idid==1 ) CYCLE
          EXIT
        END IF
        !
        IF( ndisk_com==0 ) non = numort_com + 1
        CALL DREORT(Ncomp,U(:,:,kod),V(:,kod),Yhp,Niv,W(:,non),S,P(:,non),&
          Ip(:,non),Stowa,jflag)
        !
        IF( jflag==30 ) THEN
          Iflag = 30
          !     .....................EXIT
          RETURN
          !
        ELSEIF( jflag==10 ) THEN
          xop_com = Xpts(kop_com)
          IF( ndisk_com==0 ) kod = kop_com
          !              ............EXIT
          GOTO 50
          !
        ELSEIF( jflag==0 ) THEN
          !
          !                       ************************************************
          !                           STORE ORTHONORMALIZED VECTORS INTO SOLUTION
          !                           VECTORS.
          !
          IF( numort_com>=Mxnon ) THEN
            IF( x_com/=xend_com ) THEN
              Iflag = 13
              !     .....................EXIT
              RETURN
            END IF
          END IF
          !
          numort_com = numort_com + 1
          CALL DSTOR1(Yhp(:,1),U(:,1,kod),Yhp(:,nfcp1),V(:,kod),1,ndisk_com,ntape_com)
          !
          !                       ************************************************
          !                           STORE ORTHONORMALIZATION INFORMATION,
          !                           INITIALIZE INTEGRATION FLAG, AND CONTINUE
          !                           INTEGRATION TO THE NEXT ORTHONORMALIZATION
          !                           POINT OR OUTPUT POINT.
          !
          Z(numort_com) = x_com
          IF( inhomo_com==1 .AND. nps_com==0 ) c_com = S(nfcp1)*c_com
          IF( ndisk_com/=0 ) THEN
            IF( inhomo_com==1 ) WRITE (ntape_com) (W(j,1),j=1,Nfcc)
            WRITE (ntape_com) (Ip(j,1),j=1,Nfcc), (P(j,1),j=1,Ntp)
          END IF
          info_com(1) = 0
          jon = jon + 1
          !                 ......EXIT
          IF( nopg_com==1 .AND. x_com/=xop_com ) GOTO 50
          !
          !                       ************************************************
          !                           CONTINUE INTEGRATION IF WE ARE NOT AT AN
          !                           OUTPUT POINT.
          !
          !           ............EXIT
          IF( idid/=1 ) EXIT
          !
          !                          *********************************************
          !                              CONTINUE INTEGRATION IF WE ARE NOT AT AN
          !                              OUTPUT POINT.
          !
          !           ...............EXIT
        ELSEIF( idid/=1 ) THEN
          EXIT
          !                    ......EXIT
        END IF
      ELSE
        info_com(1) = 1
        !                    ......EXIT
        IF( idid/=-1 ) THEN
          Iflag = 20 - idid
          !     .....................EXIT
          RETURN
        END IF
      END IF
    END DO
    !
    !           STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
    !           SOLUTION IN V AT THE OUTPUT POINTS.
    !
    CALL DSTOR1(U(:,1,kod),Yhp(:,1),V(:,kod),Yhp(:,nfcp1),0,ndisk_com,ntape_com)
  END DO
  !        ***************************************************************
  !        ***************************************************************
  !
  Iflag = 0
  RETURN

CONTAINS
  SUBROUTINE DBVDER_2(X,Y,Yp)
    REAL(DP), INTENT(IN) :: X
    REAL(DP), INTENT(IN) :: Y(:)
    REAL(DP), INTENT(OUT) :: Yp(:)
    !
    REAL(DP) :: g(SIZE(Y))
    !
    g = 0._SP
    CALL DBVDER(X,Y,Yp,g)
    !
  END SUBROUTINE DBVDER_2
END SUBROUTINE DRKFAB