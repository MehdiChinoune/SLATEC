!** DRKFAB
SUBROUTINE DRKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,&
    W,S,Stowa,Work,Iwork,Nfcc)
  !>
  !  Subsidiary to DBVSUP
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
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE DML, ONLY : C, INHomo, KKKint, LLLint, X, XBEg, XENd, XOP, INFo, KOP, &
    AE, RE, NOPg, NDIsk, NTApe, NEQ, INTeg, NPS, NUMort
  !
  INTEGER :: Iflag, Mxnon, Ncomp, Nfc, Nfcc, Niv, Ntp, Nxpts
  INTEGER :: Ip(Nfcc,Mxnon+1), Iwork(*)
  REAL(8) :: P(Ntp,Mxnon+1), S(Nfc+1), Stowa(:), U(Ncomp,Nfc,Nxpts), &
    V(Ncomp,Nxpts), W(Nfcc,Mxnon+1), Work(*), Xpts(Nxpts), Yhp(Ncomp,Nfc+1), Z(Mxnon+1)
  INTEGER :: idid, ipar(1), j, jflag, jon, kod, kopp, nfcp1, non
  REAL(8) :: xxop, ret(1), aet(1)
  !
  !      *****************************************************************
  !       INITIALIZATION OF COUNTERS AND VARIABLES.
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 220
  !        BEGIN BLOCK PERMITTING ...EXITS TO 10
  !* FIRST EXECUTABLE STATEMENT  DRKFAB
  kod = 1
  non = 1
  X = XBEg
  jon = 1
  INFo(1) = 0
  INFo(2) = 0
  INFo(3) = 1
  INFo(4) = 1
  Work(1) = XENd
  ipar = 0
  !        ...EXIT
  IF ( NOPg/=0 ) THEN
    INFo(3) = 0
    IF ( X==Z(1) ) jon = 2
  END IF
  nfcp1 = Nfc + 1
  !
  !        ***************************************************************
  !        *****BEGINNING OF INTEGRATION LOOP AT OUTPUT
  !        POINTS.******************
  !        ***************************************************************
  !
  DO kopp = 2, Nxpts
    KOP = kopp
    XOP = Xpts(KOP)
    IF ( NDIsk==0 ) kod = KOP
    !
    !
    !              STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
    !
    !              BEGIN BLOCK PERMITTING ...EXITS TO 190
    !                 BEGIN BLOCK PERMITTING ...EXITS TO 30
    50  xxop = XOP
    !                 ...EXIT
    IF ( NOPg/=0 ) THEN
      IF ( XENd>XBEg.AND.XOP>Z(jon) ) xxop = Z(jon)
      IF ( XENd<XBEg.AND.XOP<Z(jon) ) xxop = Z(jon)
    END IF
    DO
      !
      !                 ******************************************************
      !                    BEGIN BLOCK PERMITTING ...EXITS TO 170
      IF ( INTeg==2 ) THEN
        !                       DDEABM INTEGRATOR
        !
        ret(1) = RE
        aet(1) = AE
        CALL DDEABM(DBVDER_2,NEQ,X,Yhp,xxop,INFo,ret,aet,idid,Work,KKKint,Iwork,&
          LLLint)
      ELSE
        !                       DDERKF INTEGRATOR
        !
        ret(1) = RE
        aet(1) = AE
        CALL DDERKF(DBVDER_2,NEQ,X,Yhp,xxop,INFo,ret,aet,idid,Work,KKKint,Iwork,&
          LLLint)
      END IF
      IF ( idid>=1 ) THEN
        !
        !                       ************************************************
        !                           GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR
        !                           ORTHONORMALIZATION (TEMPORARILY USING U AND
        !                           V IN THE TEST)
        !
        IF ( NOPg==0 ) THEN
          jflag = 1
          IF ( INHomo==3.AND.X==XENd ) jflag = 3
        ELSEIF ( xxop==Z(jon) ) THEN
          jflag = 2
        ELSE
          !
          !                             ******************************************
          !                                 CONTINUE INTEGRATION IF WE ARE NOT AT
          !                                 AN OUTPUT POINT.
          !
          !           ..................EXIT
          !                    .........EXIT
          IF ( idid==1 ) CYCLE
          EXIT
        END IF
        !
        IF ( NDIsk==0 ) non = NUMort + 1
        CALL DREORT(Ncomp,U(:,:,kod),V(:,kod),Yhp,Niv,W(:,non),S,P(:,non),&
          Ip(:,non),Stowa,jflag)
        !
        IF ( jflag==30 ) THEN
          Iflag = 30
          !     .....................EXIT
          RETURN
          !
        ELSEIF ( jflag==10 ) THEN
          XOP = Xpts(KOP)
          IF ( NDIsk==0 ) kod = KOP
          !              ............EXIT
          GOTO 50
          !
        ELSEIF ( jflag==0 ) THEN
          !
          !                       ************************************************
          !                           STORE ORTHONORMALIZED VECTORS INTO SOLUTION
          !                           VECTORS.
          !
          IF ( NUMort>=Mxnon ) THEN
            IF ( X/=XENd ) THEN
              Iflag = 13
              !     .....................EXIT
              RETURN
            END IF
          END IF
          !
          NUMort = NUMort + 1
          CALL DSTOR1(Yhp(:,1),U(:,1,kod),Yhp(:,nfcp1),V(:,kod),1,NDIsk,NTApe)
          !
          !                       ************************************************
          !                           STORE ORTHONORMALIZATION INFORMATION,
          !                           INITIALIZE INTEGRATION FLAG, AND CONTINUE
          !                           INTEGRATION TO THE NEXT ORTHONORMALIZATION
          !                           POINT OR OUTPUT POINT.
          !
          Z(NUMort) = X
          IF ( INHomo==1.AND.NPS==0 ) C = S(nfcp1)*C
          IF ( NDIsk/=0 ) THEN
            IF ( INHomo==1 ) WRITE (NTApe) (W(j,1),j=1,Nfcc)
            WRITE (NTApe) (Ip(j,1),j=1,Nfcc), (P(j,1),j=1,Ntp)
          END IF
          INFo(1) = 0
          jon = jon + 1
          !                 ......EXIT
          IF ( NOPg==1.AND.X/=XOP ) GOTO 50
          !
          !                       ************************************************
          !                           CONTINUE INTEGRATION IF WE ARE NOT AT AN
          !                           OUTPUT POINT.
          !
          !           ............EXIT
          IF ( idid/=1 ) EXIT
          !
          !                          *********************************************
          !                              CONTINUE INTEGRATION IF WE ARE NOT AT AN
          !                              OUTPUT POINT.
          !
          !           ...............EXIT
        ELSEIF ( idid/=1 ) THEN
          EXIT
          !                    ......EXIT
        END IF
      ELSE
        INFo(1) = 1
        !                    ......EXIT
        IF ( idid/=-1 ) THEN
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
    CALL DSTOR1(U(:,1,kod),Yhp(:,1),V(:,kod),Yhp(:,nfcp1),0,NDIsk,NTApe)
  END DO
  !        ***************************************************************
  !        ***************************************************************
  !
  Iflag = 0
  RETURN

CONTAINS
  SUBROUTINE DBVDER_2(X,Y,Yp)
    REAL(8) :: X
    REAL(8) :: Y(:), Yp(:)
    REAL(8) :: g(SIZE(Y))

    g = 0.
    CALL DBVDER(X,Y,Yp,G)
  END SUBROUTINE DBVDER_2
END SUBROUTINE DRKFAB
