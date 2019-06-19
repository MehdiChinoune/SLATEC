!** EXBVP
SUBROUTINE EXBVP(Y,Nrowy,Xpts,A,Nrowa,Alpha,B,Nrowb,Beta,Iflag,Work,Iwork)
  !> Subsidiary to BVSUP
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
  USE ML, ONLY : ncomp_com, nfc_com, tol_com, nxpts_com, nic_com, nopg_com, &
    mxnon_com, ndisk_com, ntp_com, nfcc_com, x_com, xbeg_com, xend_com, kkkzpw_com, &
    needw_com, neediw_com, k1_com, k2_com, k3_com, k4_com, k5_com, k6_com, k7_com, &
    k8_com, k9_com, k10_com, l1_com, lpar_com
  USE service, ONLY : XERMSG
  INTEGER :: Nrowa, Nrowb, Nrowy, Iflag, Iwork(*)
  REAL(SP) :: A(:,:), Alpha(:), B(:,:), Beta(:), Work(*), Xpts(:), Y(:,:)
  INTEGER :: nsafiw, nsafw, iexp, inc, kotc
  REAL(SP) :: xl, zquit
  CHARACTER(8) :: xern1, xern2
  !* FIRST EXECUTABLE STATEMENT  EXBVP
  kotc = 1
  iexp = 0
  IF( Iwork(7)==-1 ) iexp = Iwork(8)
  DO
    !
    !     COMPUTE ORTHONORMALIZATION TOLERANCES.
    !
    tol_com = 10.0**((-lpar_com-iexp)*2)
    !
    Iwork(8) = iexp
    mxnon_com = Iwork(2)
    !
    !- *********************************************************************
    !- *********************************************************************
    !
    CALL BVPOR(Y,Nrowy,ncomp_com,Xpts,nxpts_com,A,Nrowa,Alpha,nic_com,B,Nrowb,Beta,nfc_com,&
      Iflag,Work(1),mxnon_com,Work(k1_com),ntp_com,Iwork(18),Work(k2_com),Iwork(16)&
      ,Work(k3_com),Work(k4_com),Work(k5_com),Work(k6_com),Work(k7_com),Work(k8_com:k9_com-1),&
      Work(k10_com),Iwork(l1_com),nfcc_com)
    !
    !- *********************************************************************
    !- *********************************************************************
    !     IF MGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE
    !     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE
    !     A MAXIMUM OF 2 TIMES.
    !
    IF( Iflag/=30 ) THEN
      !
      !- *********************************************************************
      !     IF BVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF
      !     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN
      !     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM
      !
      IF( Iflag==13 ) THEN
        xl = ABS(xend_com-xbeg_com)
        zquit = ABS(x_com-xbeg_com)
        inc = INT( 1.5*xl/zquit*(mxnon_com+1) )
        IF( ndisk_com/=1 ) THEN
          nsafw = inc*kkkzpw_com + needw_com
          nsafiw = inc*nfcc_com + neediw_com
        ELSE
          nsafw = needw_com + inc
          nsafiw = neediw_com
        END IF
        !
        WRITE (xern1,'(I8)') nsafw
        WRITE (xern2,'(I8)') nsafiw
        CALL XERMSG('EXBVP',&
          'IN BVSUP, PREDICTED STORAGE ALLOCATION FOR WORK ARRAY IS '&
          //xern1//', PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY IS '//xern2,1,0)
      END IF
      !
      Iwork(1) = mxnon_com
      EXIT
    ELSEIF( kotc==3 .OR. nopg_com==1 ) THEN
      Iwork(1) = mxnon_com
      EXIT
    ELSE
      kotc = kotc + 1
      iexp = iexp - 2
    END IF
  END DO
END SUBROUTINE EXBVP
