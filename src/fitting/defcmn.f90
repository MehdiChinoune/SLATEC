!** DEFCMN
SUBROUTINE DEFCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkptin,Mdein,Mdeout,&
    Coeff,Bf,Xtemp,Ptemp,Bkpt,G,Mdg,W,Mdw,Lw)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DEFC
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (EFCMN-S, DEFCMN-D)
  !***
  ! **Author:**  Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to DEFC( ).
  !     This subprogram does weighted least squares fitting of data by
  !     B-spline curves.
  !     The documentation for DEFC( ) has complete usage instructions.
  !
  !***
  ! **See also:**  DEFC
  !***
  ! **Routines called:**  DBNDAC, DBNDSL, DCOPY, DFSPVN, DSCAL, DSORT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and extensively revised (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   900604  DP version created from SP version.  (RWC)
  
  INTEGER Lw, Mdein, Mdeout, Mdg, Mdw, Nbkpt, Ndata, Nord
  REAL(8) :: Bf(Nord,*), Bkpt(*), Bkptin(*), Coeff(*), G(Mdg,*), &
    Ptemp(*), Sddata(*), W(Mdw,*), Xdata(*), Xtemp(*), &
    Ydata(*)
  !
  EXTERNAL DBNDAC, DBNDSL, DCOPY, DFSPVN, DSCAL, DSORT, XERMSG
  !
  REAL(8) :: dummy, rnorm, xmax, xmin, xval
  INTEGER i, idata, ileft, intseq, ip, ir, irow, l, mt, n, nb, &
    nordm1, nordp1, np1
  CHARACTER(8) :: xern1, xern2
  !
  !* FIRST EXECUTABLE STATEMENT  DEFCMN
  !
  !     Initialize variables and analyze input.
  !
  n = Nbkpt - Nord
  np1 = n + 1
  !
  !     Initially set all output coefficients to zero.
  !
  CALL DCOPY(n,0.D0,0,Coeff,1)
  Mdeout = -1
  IF ( Nord<1.OR.Nord>20 ) THEN
    CALL XERMSG('SLATEC','DEFCMN',&
      'IN DEFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',3,1)
    RETURN
  ENDIF
  !
  IF ( Nbkpt<2*Nord ) THEN
    CALL XERMSG('SLATEC','DEFCMN',&
      'IN DEFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.',4,1)
    RETURN
  ENDIF
  !
  IF ( Ndata<0 ) THEN
    CALL XERMSG('SLATEC','DEFCMN',&
      'IN DEFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',5,&
      1)
    RETURN
  ENDIF
  !
  nb = (Nbkpt-Nord+3)*(Nord+1) + (Nbkpt+1)*(Nord+1) + 2*MAX(Nbkpt,Ndata)&
    + Nbkpt + Nord**2
  IF ( Lw<nb ) THEN
    WRITE (xern1,'(I8)') nb
    WRITE (xern2,'(I8)') Lw
    CALL XERMSG('SLATEC','DEFCMN',&
      'IN DEFC, INSUFFICIENT STORAGE FOR W(*).  CHECK FORMULA THAT READS LW.GE. ... .  NEED = '//xern1//' GIVEN = '//&
      xern2,6,1)
    Mdeout = -1
    RETURN
  ENDIF
  !
  IF ( Mdein/=1.AND.Mdein/=2 ) THEN
    CALL XERMSG('SLATEC','DEFCMN',&
      'IN DEFC, INPUT VALUE OF MDEIN MUST BE 1-2.',7,1)
    RETURN
  ENDIF
  !
  !     Sort the breakpoints.
  !
  CALL DCOPY(Nbkpt,Bkptin,1,Bkpt,1)
  CALL DSORT(Bkpt,[dummy],Nbkpt,1)
  !
  !     Save interval containing knots.
  !
  xmin = Bkpt(Nord)
  xmax = Bkpt(np1)
  nordm1 = Nord - 1
  nordp1 = Nord + 1
  !
  !     Process least squares equations.
  !
  !     Sort data and an array of pointers.
  !
  CALL DCOPY(Ndata,Xdata,1,Xtemp,1)
  DO i = 1, Ndata
    Ptemp(i) = i
  ENDDO
  !
  IF ( Ndata>0 ) THEN
    CALL DSORT(Xtemp,Ptemp,Ndata,2)
    xmin = MIN(xmin,Xtemp(1))
    xmax = MAX(xmax,Xtemp(Ndata))
  ENDIF
  !
  !     Fix breakpoint array if needed. This should only involve very
  !     minor differences with the input array of breakpoints.
  !
  DO i = 1, Nord
    Bkpt(i) = MIN(Bkpt(i),xmin)
  ENDDO
  !
  DO i = np1, Nbkpt
    Bkpt(i) = MAX(Bkpt(i),xmax)
  ENDDO
  !
  !     Initialize parameters of banded matrix processor, DBNDAC( ).
  !
  mt = 0
  ip = 1
  ir = 1
  ileft = Nord
  intseq = 1
  DO idata = 1, Ndata
    !
    !        Sorted indices are in PTEMP(*).
    !
    l = INT( Ptemp(idata) )
    xval = Xdata(l)
    !
    !        When interval changes, process equations in the last block.
    !
    IF ( xval>=Bkpt(ileft+1) ) THEN
      CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
      mt = 0
      !
      !           Move pointer up to have BKPT(ILEFT).LE.XVAL, ILEFT.LE.N.
      !
      DO ileft = ileft, n
        IF ( xval<Bkpt(ileft+1) ) EXIT
        IF ( Mdein==2 ) THEN
          !
          !                 Data is being sequentially accumulated.
          !                 Transfer previously accumulated rows from W(*,*) to
          !                 G(*,*) and process them.
          !
          CALL DCOPY(nordp1,W(intseq,1),Mdw,G(ir,1),Mdg)
          CALL DBNDAC(G,Mdg,Nord,ip,ir,1,intseq)
          intseq = intseq + 1
        ENDIF
      ENDDO
    ENDIF
    !
    !        Obtain B-spline function value.
    !
    CALL DFSPVN(Bkpt,Nord,1,xval,ileft,Bf)
    !
    !        Move row into place.
    !
    irow = ir + mt
    mt = mt + 1
    CALL DCOPY(Nord,Bf,1,G(irow,1),Mdg)
    G(irow,nordp1) = Ydata(l)
    !
    !        Scale data if uncertainty is nonzero.
    !
    IF ( Sddata(l)/=0.D0 ) CALL DSCAL(nordp1,1.D0/Sddata(l),G(irow,1),Mdg)
    !
    !        When staging work area is exhausted, process rows.
    !
    IF ( irow==Mdg-1 ) THEN
      CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
      mt = 0
    ENDIF
  ENDDO
  !
  !     Process last block of equations.
  !
  CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
  !
  !     Finish processing any previously accumulated rows from W(*,*)
  !     to G(*,*).
  !
  IF ( Mdein==2 ) THEN
    DO i = intseq, np1
      CALL DCOPY(nordp1,W(i,1),Mdw,G(ir,1),Mdg)
      CALL DBNDAC(G,Mdg,Nord,ip,ir,1,MIN(n,i))
    ENDDO
  ENDIF
  !
  !     Last call to adjust block positioning.
  !
  CALL DCOPY(nordp1,0.D0,0,G(ir,1),Mdg)
  CALL DBNDAC(G,Mdg,Nord,ip,ir,1,np1)
  !
  !     Transfer accumulated rows from G(*,*) to W(*,*) for
  !     possible later sequential accumulation.
  !
  DO i = 1, np1
    CALL DCOPY(nordp1,G(i,1),Mdg,W(i,1),Mdw)
  ENDDO
  !
  !     Solve for coefficients when possible.
  !
  DO i = 1, n
    IF ( G(i,1)==0.D0 ) THEN
      Mdeout = 2
      RETURN
    ENDIF
  ENDDO
  !
  !     All the diagonal terms in the accumulated triangular
  !     matrix are nonzero.  The solution can be computed but
  !     it may be unsuitable for further use due to poor
  !     conditioning or the lack of constraints.  No checking
  !     for either of these is done here.
  !
  CALL DBNDSL(1,G,Mdg,Nord,ip,ir,Coeff,n,rnorm)
  Mdeout = 1
END SUBROUTINE DEFCMN
