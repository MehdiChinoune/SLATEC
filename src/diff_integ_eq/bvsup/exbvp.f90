!** EXBVP
SUBROUTINE EXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  USE ML, ONLY : NCOmp, NFC, TOL, NXPts, NIC, NOPg, MXNon, NDIsk, NTP, NFCc, X, &
    XBEg, XENd, KKKzpw, NEEdw, NEEdiw, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, L1, LPAr
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

  INTEGER Nrowa, Nrowb, Nrowy, nsafiw, nsafw, iexp, Iflag, inc, Iwork(*), kotc
  REAL A(Nrowa,*), Alpha(*), B(Nrowb,*), Beta(*), Work(*), xl, Xpts(*), Y(Nrowy,*), zquit
  CHARACTER(8) :: xern1, xern2
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
        END IF
        !
        WRITE (xern1,'(I8)') nsafw
        WRITE (xern2,'(I8)') nsafiw
        CALL XERMSG('SLATEC','EXBVP',&
          'IN BVSUP, PREDICTED STORAGE ALLOCATION FOR WORK ARRAY IS '&
          //xern1//', PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY IS '//xern2,1,0)
      END IF
      !
      Iwork(1) = MXNon
      EXIT
    ELSEIF ( kotc==3.OR.NOPg==1 ) THEN
      Iwork(1) = MXNon
      EXIT
    ELSE
      kotc = kotc + 1
      iexp = iexp - 2
    END IF
  END DO
END SUBROUTINE EXBVP
