!DECK RKFAB
SUBROUTINE RKFAB(Ncomp,Xpts,Nxpts,Nfc,Iflag,Z,Mxnon,P,Ntp,Ip,Yhp,Niv,U,V,&
    W,S,Stowa,G,Work,Iwork,Nfcc)
  IMPLICIT NONE
  REAL AE, C, G, P, PWCnd, PX, RE, S, Stowa, TND, TOL, &
    U, V, W, Work, X, XBEg, XENd, XOP
  REAL XOT, Xpts, XSAv, xxop, Yhp, Z
  INTEGER ICOco, idid, Iflag, IGOfx, INDpvt, INFo, INHomo, INTeg, &
    Ip, ipar, ISTkop, IVP, Iwork, j, jflag, jon, K1, K10, &
    K11, K2
  INTEGER K3, K4, K5, K6, K7, K8, K9, KKKint, KKKzpw, KNSwot, &
    kod, KOP, kopp, L1, L2, LLLint, LOTjp, MNSwot, Mxnon, &
    MXNond
  INTEGER Ncomp, NCOmpd, NDIsk, NEEdiw, NEEdw, NEQ, NEQivp, Nfc, &
    Nfcc, NFCcd, NFCd, nfcp1, NIC, Niv, non, NOPg, NPS, &
    NSWot, NTApe, Ntp
  INTEGER NTPd, NUMort, Nxpts, NXPtsd
  !***BEGIN PROLOGUE  RKFAB
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (RKFAB-S, DRKFAB-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  ! **********************************************************************
  !
  !     Subroutine RKFAB integrates the initial value equations using
  !     the variable-step RUNGE-KUTTA-FEHLBERG integration scheme or
  !     the variable-order ADAMS method and orthonormalization
  !     determined by a linear dependence test.
  !
  ! **********************************************************************
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  BVDER, DEABM, DERKF, REORT, STOR1
  !***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML8SZ
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  RKFAB
  !
  DIMENSION P(Ntp,*), Ip(Nfcc,*), U(Ncomp,Nfc,*), V(Ncomp,*), W(Nfcc,*)&
    , Z(*), Yhp(Ncomp,*), Xpts(*), S(*), Stowa(*), Work(*), &
    Iwork(*), G(*)
  !
  ! **********************************************************************
  !
  COMMON /ML8SZ / C, XSAv, IGOfx, INHomo, IVP, NCOmpd, NFCd
  COMMON /ML15TO/ PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, INFo(15)&
    , ISTkop, KNSwot, KOP, LOTjp, MNSwot, NSWot
  COMMON /ML18JR/ AE, RE, TOL, NXPtsd, NIC, NOPg, MXNond, NDIsk, &
    NTApe, NEQ, INDpvt, INTeg, NPS, NTPd, NEQivp, &
    NUMort, NFCcd, ICOco
  COMMON /ML17BW/ KKKzpw, NEEdw, NEEdiw, K1, K2, K3, K4, K5, K6, &
    K7, K8, K9, K10, K11, L1, L2, KKKint, LLLint
  !
  EXTERNAL BVDER
  !
  ! **********************************************************************
  !  INITIALIZATION OF COUNTERS AND VARIABLES.
  !
  !***FIRST EXECUTABLE STATEMENT  RKFAB
  kod = 1
  non = 1
  X = XBEg
  jon = 1
  INFo(1) = 0
  INFo(2) = 0
  INFo(3) = 1
  INFo(4) = 1
  Work(1) = XENd
  IF ( NOPg/=0 ) THEN
    INFo(3) = 0
    IF ( X==Z(1) ) jon = 2
  ENDIF
  nfcp1 = Nfc + 1
  !
  ! **********************************************************************
  ! *****BEGINNING OF INTEGRATION LOOP AT OUTPUT POINTS.******************
  ! **********************************************************************
  !
  DO kopp = 2, Nxpts
    KOP = kopp
    !
    50     XOP = Xpts(KOP)
    IF ( NDIsk==0 ) kod = KOP
    !
    !     STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
    !
    100    xxop = XOP
    IF ( NOPg/=0 ) THEN
      IF ( XENd>XBEg.AND.XOP>Z(jon) ) xxop = Z(jon)
      IF ( XENd<XBEg.AND.XOP<Z(jon) ) xxop = Z(jon)
    ENDIF
    !
    ! **********************************************************************
    150    IF ( INTeg==2 ) THEN
    !     DEABM INTEGRATOR
    !
    CALL DEABM(BVDER,NEQ,X,Yhp,xxop,INFo,[RE],[AE],idid,Work,KKKint,Iwork,&
      LLLint,G,[ipar])
  ELSE
    !     DERKF INTEGRATOR
    !
    CALL DERKF(BVDER,NEQ,X,Yhp,xxop,INFo,[RE],[AE],idid,Work,KKKint,Iwork,&
      LLLint,G,[ipar])
  ENDIF
  IF ( idid>=1 ) THEN
    !
    ! **********************************************************************
    !     GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR ORTHONORMALIZATION
    !     (TEMPORARILY USING U AND V IN THE TEST)
    !
    IF ( NOPg==0 ) THEN
      jflag = 1
      IF ( INHomo==3.AND.X==XENd ) jflag = 3
    ELSE
      IF ( xxop/=Z(jon) ) GOTO 200
      jflag = 2
    ENDIF
    !
    IF ( NDIsk==0 ) non = NUMort + 1
    CALL REORT(Ncomp,U(1,1,kod),V(1,kod),Yhp,Niv,W(1,non),S,P(1,non),&
      Ip(1,non),Stowa,jflag)
    !
    IF ( jflag/=30 ) THEN
      !
      IF ( jflag==10 ) GOTO 50
      !
      IF ( jflag==0 ) THEN
        !
        ! **********************************************************************
        !     STORE ORTHONORMALIZED VECTORS INTO SOLUTION VECTORS.
        !
        IF ( NUMort>=Mxnon ) THEN
          IF ( X/=XENd ) THEN
            Iflag = 13
            RETURN
          ENDIF
        ENDIF
        !
        NUMort = NUMort + 1
        CALL STOR1(Yhp,U(1,1,kod),Yhp(1,nfcp1),V(1,kod),1,NDIsk,NTApe)
        !
        ! **********************************************************************
        !     STORE ORTHONORMALIZATION INFORMATION, INITIALIZE
        !     INTEGRATION FLAG, AND CONTINUE INTEGRATION TO THE NEXT
        !     ORTHONORMALIZATION POINT OR OUTPUT POINT.
        !
        Z(NUMort) = X
        IF ( INHomo==1.AND.NPS==0 ) C = S(nfcp1)*C
        IF ( NDIsk/=0 ) THEN
          IF ( INHomo==1 ) WRITE (NTApe) (W(j,1),j=1,Nfcc)
          WRITE (NTApe) (Ip(j,1),j=1,Nfcc), (P(j,1),j=1,Ntp)
        ENDIF
        INFo(1) = 0
        jon = jon + 1
        IF ( NOPg==1.AND.X/=XOP ) GOTO 100
      ENDIF
    ELSE
      Iflag = 30
      RETURN
    ENDIF
  ELSE
    INFo(1) = 1
    IF ( idid==-1 ) GOTO 150
    Iflag = 20 - idid
    RETURN
  ENDIF
  !
  ! **********************************************************************
  !     CONTINUE INTEGRATION IF WE ARE NOT AT AN OUTPUT POINT.
  !
  200    IF ( idid==1 ) GOTO 150
  !
  !     STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
  !     SOLUTION IN V AT THE OUTPUT POINTS.
  !
  CALL STOR1(U(1,1,kod),Yhp,V(1,kod),Yhp(1,nfcp1),0,NDIsk,NTApe)
  ENDDO
  ! **********************************************************************
  ! **********************************************************************
  !
  Iflag = 0
END SUBROUTINE RKFAB
