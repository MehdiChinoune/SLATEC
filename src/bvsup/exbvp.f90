!** EXBVP
SUBROUTINE EXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (EXBVP-S, DEXBVP-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine is used to execute the basic technique for solving
  !  the two-point boundary value problem
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  BVPOR, XERMSG
  !***
  ! COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890921  Realigned order of variables in certain COMMON blocks.
  !           (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   910722  Updated AUTHOR section.  (ALS)

  INTEGER NEQivp, NFC, NFCc, NIC, NOPg, NPS, Nrowa, Nrowb, Nrowy, &
    nsafiw, nsafw, NSWot, NTApe, NTP, NUMort, NXPts
  REAL A(Nrowa,*), AE, Alpha(*), B(Nrowb,*), Beta(*), C, EPS, FOUru, PWCnd, PX, RE, &
    SQOvfl, SRU, TND, TOL, TWOu, URO, Work(*), X, XBEg
  REAL XENd, xl, XOP, XOT, Xpts(*), XSAv, Y(Nrowy,*), zquit
  INTEGER ICOco, iexp, Iflag, IGOfx, inc, INDpvt, INFo, INHomo, &
    INTeg, ISTkop, IVP, Iwork(*), K1, K10, K11, K2, K3, K4, K5, K6
  INTEGER K7, K8, K9, KKKint, KKKzpw, KNSwot, KOP, kotc, L1, L2, &
    LLLint, LOTjp, LPAr, MNSwot, MXNon, NCOmp, NDIsk, NEEdiw, NEEdw, NEQ
  CHARACTER(8) :: xern1, xern2
  !
  !     ****************************************************************
  !
  COMMON /ML8SZ / C, XSAv, IGOfx, INHomo, IVP, NCOmp, NFC
  COMMON /ML18JR/ AE, RE, TOL, NXPts, NIC, NOPg, MXNon, NDIsk, &
    NTApe, NEQ, INDpvt, INTeg, NPS, NTP, NEQivp, NUMort, NFCc, ICOco
  COMMON /ML15TO/ PX, PWCnd, TND, X, XBEg, XENd, XOT, XOP, INFo(15), &
    ISTkop, KNSwot, KOP, LOTjp, MNSwot, NSWot
  COMMON /ML17BW/ KKKzpw, NEEdw, NEEdiw, K1, K2, K3, K4, K5, K6, &
    K7, K8, K9, K10, K11, L1, L2, KKKint, LLLint
  !
  COMMON /ML5MCO/ URO, SRU, EPS, SQOvfl, TWOu, FOUru, LPAr
  !
  !* FIRST EXECUTABLE STATEMENT  EXBVP
  kotc = 1
  iexp = 0
  IF ( Iwork(7)==-1 ) iexp = Iwork(8)
  DO
    !
    !     COMPUTE ORTHONORMALIZATION TOLERANCES.
    !
    TOL = 10.0**((-LPAr-iexp)*2)
    !
    Iwork(8) = iexp
    MXNon = Iwork(2)
    !
    !- *********************************************************************
    !- *********************************************************************
    !
    CALL BVPOR(Y,Nrowy,NCOmp,Xpts,NXPts,A,Nrowa,Alpha,NIC,B,Nrowb,Beta,NFC,&
      Iflag,Work(1),MXNon,Work(K1),NTP,Iwork(18),Work(K2),Iwork(16)&
      ,Work(K3),Work(K4),Work(K5),Work(K6),Work(K7),Work(K8),&
      Work(K9),Work(K10),Iwork(L1),NFCc)
    !
    !- *********************************************************************
    !- *********************************************************************
    !     IF MGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE
    !     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE
    !     A MAXIMUM OF 2 TIMES.
    !
    IF ( Iflag/=30 ) THEN
      !
      !- *********************************************************************
      !     IF BVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF
      !     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN
      !     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM
      !
      IF ( Iflag==13 ) THEN
        xl = ABS(XENd-XBEg)
        zquit = ABS(X-XBEg)
        inc = INT( 1.5*xl/zquit*(MXNon+1) )
        IF ( NDIsk/=1 ) THEN
          nsafw = inc*KKKzpw + NEEdw
          nsafiw = inc*NFCc + NEEdiw
        ELSE
          nsafw = NEEdw + inc
          nsafiw = NEEdiw
        ENDIF
        !
        WRITE (xern1,'(I8)') nsafw
        WRITE (xern2,'(I8)') nsafiw
        CALL XERMSG('SLATEC','EXBVP',&
          'IN BVSUP, PREDICTED STORAGE ALLOCATION FOR WORK ARRAY IS '&
          //xern1//', PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY IS '//xern2,1,0)
      ENDIF
      !
      Iwork(1) = MXNon
      EXIT
    ELSEIF ( kotc==3.OR.NOPg==1 ) THEN
      Iwork(1) = MXNon
      EXIT
    ELSE
      kotc = kotc + 1
      iexp = iexp - 2
    ENDIF
  ENDDO
END SUBROUTINE EXBVP
