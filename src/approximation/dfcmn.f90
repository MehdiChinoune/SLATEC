!** DFCMN
SUBROUTINE DFCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkptin,Nconst,Xconst,&
    Yconst,Nderiv,Mode,Coeff,Bf,Xtemp,Ptemp,Bkpt,G,Mdg,W,Mdw,Work,Iwork)
  !>
  !***
  !  Subsidiary to FC
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (FCMN-S, DFCMN-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to DFC( ).
  !     The documentation for DFC( ) has complete usage instructions.
  !
  !***
  ! **See also:**  DFC
  !***
  ! **Routines called:**  DAXPY, DBNDAC, DBNDSL, DCOPY, DFSPVD, DFSPVN,
  !                    DLSEI, DSCAL, DSORT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and extensively revised (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   900604  DP version created from SP version.  (RWC)

  INTEGER Iwork(*), Mdg, Mdw, Mode, Nbkpt, Nconst, Ndata, Nderiv(*), Nord
  REAL(8) :: Bf(Nord,*), Bkpt(*), Bkptin(*), Coeff(*), G(Mdg,*), &
    Ptemp(*), Sddata(*), W(Mdw,*), Work(*), Xconst(*), &
    Xdata(*), Xtemp(*), Yconst(*), Ydata(*)
  !
  REAL(8) :: dummy(1), prgopt(10), rnorm, rnorme, rnorml, xmax, xmin, xval, yval
  INTEGER i, idata, ideriv, ileft, intrvl, intw1, ip, ir, irow, &
    itype, iw1, iw2, l, lw, mt, n, nb, neqcon, nincon, nordm1, nordp1, np1
  LOGICAL band, new, var
  CHARACTER(8) :: xern1
  !
  !* FIRST EXECUTABLE STATEMENT  DFCMN
  !
  !     Analyze input.
  !
  dummy = 0.D0
  IF ( Nord<1.OR.Nord>20 ) THEN
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',2,1)
    Mode = -1
    RETURN
    !
  ELSEIF ( Nbkpt<2*Nord ) THEN
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.',2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF ( Ndata<0 ) THEN
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',2,1)
    Mode = -1
    RETURN
  END IF
  !
  !     Amount of storage allocated for W(*), IW(*).
  !
  iw1 = Iwork(1)
  iw2 = Iwork(2)
  nb = (Nbkpt-Nord+3)*(Nord+1) + 2*MAX(Ndata,Nbkpt) + Nbkpt + Nord**2
  !
  !     See if sufficient storage has been allocated.
  !
  IF ( iw1<nb ) THEN
    WRITE (xern1,'(I8)') nb
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = '//xern1,2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF ( Mode==1 ) THEN
    band = .TRUE.
    var = .FALSE.
    new = .TRUE.
  ELSEIF ( Mode==2 ) THEN
    band = .FALSE.
    var = .TRUE.
    new = .TRUE.
  ELSEIF ( Mode==3 ) THEN
    band = .TRUE.
    var = .FALSE.
    new = .FALSE.
  ELSEIF ( Mode==4 ) THEN
    band = .FALSE.
    var = .TRUE.
    new = .FALSE.
  ELSE
    CALL XERMSG('SLATEC','DFCMN','IN DFC, INPUT VALUE OF MODE MUST BE 1-4.',2,1)
    Mode = -1
    RETURN
  END IF
  Mode = 0
  !
  !     Sort the breakpoints.
  !
  CALL DCOPY(Nbkpt,Bkptin,1,Bkpt,1)
  CALL DSORT(Bkpt,dummy,Nbkpt,1)
  !
  !     Initialize variables.
  !
  neqcon = 0
  nincon = 0
  DO i = 1, Nconst
    l = Nderiv(i)
    itype = MOD(l,4)
    IF ( itype<2 ) THEN
      nincon = nincon + 1
    ELSE
      neqcon = neqcon + 1
    END IF
  END DO
  !
  !     Compute the number of variables.
  !
  n = Nbkpt - Nord
  np1 = n + 1
  lw = nb + (np1+Nconst)*np1 + 2*(neqcon+np1) + (nincon+np1) + (nincon+2)*(np1+6)
  intw1 = nincon + 2*np1
  !
  !     Save interval containing knots.
  !
  xmin = Bkpt(Nord)
  xmax = Bkpt(np1)
  !
  !     Find the smallest referenced independent variable value in any
  !     constraint.
  !
  DO i = 1, Nconst
    xmin = MIN(xmin,Xconst(i))
    xmax = MAX(xmax,Xconst(i))
  END DO
  nordm1 = Nord - 1
  nordp1 = Nord + 1
  !
  !     Define the option vector PRGOPT(1-10) for use in DLSEI( ).
  !
  prgopt(1) = 4
  !
  !     Set the covariance matrix computation flag.
  !
  prgopt(2) = 1
  IF ( var ) THEN
    prgopt(3) = 1
  ELSE
    prgopt(3) = 0
  END IF
  !
  !     Increase the rank determination tolerances for both equality
  !     constraint equations and least squares equations.
  !
  prgopt(4) = 7
  prgopt(5) = 4
  prgopt(6) = 1.D-4
  !
  prgopt(7) = 10
  prgopt(8) = 5
  prgopt(9) = 1.D-4
  !
  prgopt(10) = 1
  !
  !     Turn off work array length checking in DLSEI( ).
  !
  Iwork(1) = 0
  Iwork(2) = 0
  !
  !     Initialize variables and analyze input.
  !
  IF ( new ) THEN
    !
    !        To process least squares equations sort data and an array of
    !        pointers.
    !
    CALL DCOPY(Ndata,Xdata,1,Xtemp,1)
    DO i = 1, Ndata
      Ptemp(i) = i
    END DO
    !
    IF ( Ndata>0 ) THEN
      CALL DSORT(Xtemp,Ptemp,Ndata,2)
      xmin = MIN(xmin,Xtemp(1))
      xmax = MAX(xmax,Xtemp(Ndata))
    END IF
    !
    !        Fix breakpoint array if needed.
    !
    DO i = 1, Nord
      Bkpt(i) = MIN(Bkpt(i),xmin)
    END DO
    !
    DO i = np1, Nbkpt
      Bkpt(i) = MAX(Bkpt(i),xmax)
    END DO
    !
    !        Initialize parameters of banded matrix processor, DBNDAC( ).
    !
    mt = 0
    ip = 1
    ir = 1
    ileft = Nord
    DO idata = 1, Ndata
      !
      !           Sorted indices are in PTEMP(*).
      !
      l = INT( Ptemp(idata) )
      xval = Xdata(l)
      !
      !           When interval changes, process equations in the last block.
      !
      IF ( xval>=Bkpt(ileft+1) ) THEN
        CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
        !
        !              Move pointer up to have BKPT(ILEFT).LE.XVAL,
        !                 ILEFT.LT.NP1.
        !
        DO WHILE ( xval>=Bkpt(ileft+1).AND.ileft<n )
          ileft = ileft + 1
        END DO
      END IF
      !
      !           Obtain B-spline function value.
      !
      CALL DFSPVN(Bkpt,Nord,1,xval,ileft,Bf)
      !
      !           Move row into place.
      !
      irow = ir + mt
      mt = mt + 1
      CALL DCOPY(Nord,Bf,1,G(irow,1),Mdg)
      G(irow,nordp1) = Ydata(l)
      !
      !           Scale data if uncertainty is nonzero.
      !
      IF ( Sddata(l)/=0.D0 ) CALL DSCAL(nordp1,1.D0/Sddata(l),G(irow,1),Mdg)
      !
      !           When staging work area is exhausted, process rows.
      !
      IF ( irow==Mdg-1 ) THEN
        CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
      END IF
    END DO
    !
    !        Process last block of equations.
    !
    CALL DBNDAC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
    !
    !        Last call to adjust block positioning.
    !
    G(ir,1:nordp1) = 0.D0
    CALL DBNDAC(G,Mdg,Nord,ip,ir,1,np1)
  END IF
  !
  band = band .AND. Nconst==0
  DO i = 1, n
    band = band .AND. G(i,1)/=0.D0
  END DO
  !
  !     Process banded least squares equations.
  !
  IF ( band ) THEN
    CALL DBNDSL(1,G,Mdg,Nord,ip,ir,Coeff,n,rnorm)
    RETURN
  END IF
  !
  !     Check further for sufficient storage in working arrays.
  !
  IF ( iw1<lw ) THEN
    WRITE (xern1,'(I8)') lw
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = '//xern1,2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF ( iw2<intw1 ) THEN
    WRITE (xern1,'(I8)') intw1
    CALL XERMSG('SLATEC','DFCMN',&
      'IN DFC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = '//xern1,2,1)
    Mode = -1
    RETURN
  END IF
  !
  !     Write equality constraints.
  !     Analyze constraint indicators for an equality constraint.
  !
  neqcon = 0
  DO idata = 1, Nconst
    l = Nderiv(idata)
    itype = MOD(l,4)
    IF ( itype>1 ) THEN
      ideriv = l/4
      neqcon = neqcon + 1
      ileft = Nord
      xval = Xconst(idata)
      !
      DO WHILE ( xval>=Bkpt(ileft+1).AND.ileft<n )
        ileft = ileft + 1
      END DO
      !
      CALL DFSPVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      W(neqcon,1:np1) = 0.D0
      CALL DCOPY(Nord,Bf(1,ideriv+1),1,W(neqcon,ileft-nordm1),Mdw)
      !
      IF ( itype==2 ) THEN
        W(neqcon,np1) = Yconst(idata)
      ELSE
        ileft = Nord
        yval = Yconst(idata)
        !
        DO WHILE ( yval>=Bkpt(ileft+1).AND.ileft<n )
          ileft = ileft + 1
        END DO
        !
        CALL DFSPVD(Bkpt,Nord,yval,ileft,Bf,ideriv+1)
        CALL DAXPY(Nord,-1.D0,Bf(1,ideriv+1),1,W(neqcon,ileft-nordm1),Mdw)
      END IF
    END IF
  END DO
  !
  !     Transfer least squares data.
  !
  DO i = 1, np1
    irow = i + neqcon
    W(irow,1:n) = 0.D0
    CALL DCOPY(MIN(np1-i,Nord),G(i,1),Mdg,W(irow,i),Mdw)
    W(irow,np1) = G(i,nordp1)
  END DO
  !
  !     Write inequality constraints.
  !     Analyze constraint indicators for inequality constraints.
  !
  nincon = 0
  DO idata = 1, Nconst
    l = Nderiv(idata)
    itype = MOD(l,4)
    IF ( itype<2 ) THEN
      ideriv = l/4
      nincon = nincon + 1
      ileft = Nord
      xval = Xconst(idata)
      !
      DO WHILE ( xval>=Bkpt(ileft+1).AND.ileft<n )
        ileft = ileft + 1
      END DO
      !
      CALL DFSPVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      irow = neqcon + np1 + nincon
      W(irow,1:n) = 0.D0
      intrvl = ileft - nordm1
      CALL DCOPY(Nord,Bf(1,ideriv+1),1,W(irow,intrvl),Mdw)
      !
      IF ( itype==1 ) THEN
        W(irow,np1) = Yconst(idata)
      ELSE
        W(irow,np1) = -Yconst(idata)
        CALL DSCAL(Nord,-1.D0,W(irow,intrvl),Mdw)
      END IF
    END IF
  END DO
  !
  !     Solve constrained least squares equations.
  !
  CALL DLSEI(W,Mdw,neqcon,np1,nincon,n,prgopt,Coeff,rnorme,rnorml,Mode,Work,&
    Iwork)
END SUBROUTINE DFCMN
