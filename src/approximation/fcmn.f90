!** FCMN
SUBROUTINE FCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkptin,Nconst,Xconst,&
    Yconst,Nderiv,Mode,Coeff,Bf,Xtemp,Ptemp,Bkpt,G,Mdg,W,Mdw,Work,Iwork)
  !> Subsidiary to FC
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (FCMN-S, DFCMN-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This is a companion subprogram to FC( ).
  !     The documentation for FC( ) has complete usage instructions.
  !
  !***
  ! **See also:**  FC
  !***
  ! **Routines called:**  BNDACC, BNDSOL, BSPLVD, BSPLVN, LSEI, SAXPY, SCOPY,
  !                    SSCAL, SSORT, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890618  Completely restructured and extensively revised (WRB & RWC)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900328  Added TYPE section.  (WRB)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  USE service, ONLY : XERMSG
  USE blas, ONLY : SAXPY
  USE linear, ONLY : BNDSOL, BNDACC
  USE data_handling, ONLY : SSORT
  INTEGER :: Mdg, Mdw, Mode, Nbkpt, Nconst, Ndata, Nord, Iwork(:), Nderiv(Nconst)
  REAL(SP) :: Bf(Nord,Nord), Bkpt(Nbkpt), Bkptin(:), Coeff(:), G(Mdg,Nord+1), &
    Ptemp(MAX(Nbkpt,Ndata)), Sddata(Ndata), W(Mdw,Nbkpt-Nord+1), Work(*), &
    Xconst(Nconst), Xdata(Ndata), Xtemp(MAX(Nbkpt,Ndata)), Yconst(Nconst), Ydata(Ndata)
  REAL(SP) :: dummy(1), prgopt(10), rnorm, rnorme, rnorml, xmax, xmin, xval, yval
  INTEGER :: i, idata, ideriv, ileft, intrvl, intw1, ip, ir, irow, itype, iw1, iw2, &
    l, lw, mt, n, nb, neqcon, nincon, nordm1, nordp1, np1
  LOGICAL :: band, new, var
  CHARACTER(8) :: xern1
  !
  !* FIRST EXECUTABLE STATEMENT  FCMN
  !
  !     Analyze input.
  !
  dummy = 0.
  IF( Nord<1 .OR. Nord>20 ) THEN
    CALL XERMSG('FCMN',&
      'IN FC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',2,1)
    Mode = -1
    RETURN
    !
  ELSEIF( Nbkpt<2*Nord ) THEN
    CALL XERMSG('FCMN',&
      'IN FC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.',2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF( Ndata<0 ) THEN
    CALL XERMSG('FCMN',&
      'IN FC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',2,1)
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
  IF( iw1<nb ) THEN
    WRITE (xern1,'(I8)') nb
    CALL XERMSG('FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = '//xern1,2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF( Mode==1 ) THEN
    band = .TRUE.
    var = .FALSE.
    new = .TRUE.
  ELSEIF( Mode==2 ) THEN
    band = .FALSE.
    var = .TRUE.
    new = .TRUE.
  ELSEIF( Mode==3 ) THEN
    band = .TRUE.
    var = .FALSE.
    new = .FALSE.
  ELSEIF( Mode==4 ) THEN
    band = .FALSE.
    var = .TRUE.
    new = .FALSE.
  ELSE
    CALL XERMSG('FCMN','IN FC, INPUT VALUE OF MODE MUST BE 1-4.',2,1)
    Mode = -1
    RETURN
  END IF
  Mode = 0
  !
  !     Sort the breakpoints.
  !
  Bkpt(1:Nbkpt) = Bkptin(1:Nbkpt)
  CALL SSORT(Bkpt,dummy,Nbkpt,1)
  !
  !     Initialize variables.
  !
  neqcon = 0
  nincon = 0
  DO i = 1, Nconst
    l = Nderiv(i)
    itype = MOD(l,4)
    IF( itype<2 ) THEN
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
  !     Define the option vector PRGOPT(1-10) for use in LSEI( ).
  !
  prgopt(1) = 4
  !
  !     Set the covariance matrix computation flag.
  !
  prgopt(2) = 1
  IF( var ) THEN
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
  prgopt(6) = 1.E-4
  !
  prgopt(7) = 10
  prgopt(8) = 5
  prgopt(9) = 1.E-4
  !
  prgopt(10) = 1
  !
  !     Turn off work array length checking in LSEI( ).
  !
  Iwork(1) = 0
  Iwork(2) = 0
  !
  !     Initialize variables and analyze input.
  !
  IF( new ) THEN
    !
    !        To process least squares equations sort data and an array of
    !        pointers.
    !
    Xtemp(1:Ndata) = Xdata(1:Ndata)
    DO i = 1, Ndata
      Ptemp(i) = i
    END DO
    !
    IF( Ndata>0 ) THEN
      CALL SSORT(Xtemp,Ptemp,Ndata,2)
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
    !        Initialize parameters of banded matrix processor, BNDACC( ).
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
      IF( xval>=Bkpt(ileft+1) ) THEN
        CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
        !
        !              Move pointer up to have BKPT(ILEFT)<=XVAL,
        !                 ILEFT<NP1.
        !
        DO WHILE( xval>=Bkpt(ileft+1) .AND. ileft<n )
          ileft = ileft + 1
        END DO
      END IF
      !
      !           Obtain B-spline function value.
      !
      CALL BSPLVN(Bkpt,Nord,1,xval,ileft,Bf)
      !
      !           Move row into place.
      !
      irow = ir + mt
      mt = mt + 1
      G(irow,1:Nord) = Bf(1:Nord,1)
      G(irow,nordp1) = Ydata(l)
      !
      !           Scale data if uncertainty is nonzero.
      !
      IF( Sddata(l)/=0.E0 ) G(irow,1:nordp1) = G(irow,1:nordp1)/Sddata(l)
      !
      !           When staging work area is exhausted, process rows.
      !
      IF( irow==Mdg-1 ) THEN
        CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
      END IF
    END DO
    !
    !        Process last block of equations.
    !
    CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
    !
    !        Last call to adjust block positioning.
    !
    G(ir,1:nordp1) = 0.E0
    CALL BNDACC(G,Mdg,Nord,ip,ir,1,np1)
  END IF
  !
  band = band .AND. Nconst==0
  DO i = 1, n
    band = band .AND. G(i,1)/=0.E0
  END DO
  !
  !     Process banded least squares equations.
  !
  IF( band ) THEN
    CALL BNDSOL(1,G,Mdg,Nord,ip,ir,Coeff,n,rnorm)
    RETURN
  END IF
  !
  !     Check further for sufficient storage in working arrays.
  !
  IF( iw1<lw ) THEN
    WRITE (xern1,'(I8)') lw
    CALL XERMSG('FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = '//xern1,2,1)
    Mode = -1
    RETURN
  END IF
  !
  IF( iw2<intw1 ) THEN
    WRITE (xern1,'(I8)') intw1
    CALL XERMSG('FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = '//xern1,2,1)
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
    IF( itype>1 ) THEN
      ideriv = l/4
      neqcon = neqcon + 1
      ileft = Nord
      xval = Xconst(idata)
      !
      DO WHILE( xval>=Bkpt(ileft+1) .AND. ileft<n )
        ileft = ileft + 1
      END DO
      !
      CALL BSPLVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      W(neqcon,1:np1) = 0.E0
      W(neqcon,ileft-nordm1:ileft) = Bf(1:Nord,ideriv+1)
      !
      IF( itype==2 ) THEN
        W(neqcon,np1) = Yconst(idata)
      ELSE
        ileft = Nord
        yval = Yconst(idata)
        !
        DO WHILE( yval>=Bkpt(ileft+1) .AND. ileft<n )
          ileft = ileft + 1
        END DO
        !
        CALL BSPLVD(Bkpt,Nord,yval,ileft,Bf,ideriv+1)
        CALL SAXPY(Nord,-1.E0,Bf(1,ideriv+1),1,W(neqcon,ileft-nordm1),Mdw)
      END IF
    END IF
  END DO
  !
  !     Transfer least squares data.
  !
  DO i = 1, np1
    irow = i + neqcon
    W(irow,1:n) = 0.E0
    W(irow,i:i+MIN(np1-i,Nord)-1) = G(i,1:MIN(np1-i,Nord))
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
    IF( itype<2 ) THEN
      ideriv = l/4
      nincon = nincon + 1
      ileft = Nord
      xval = Xconst(idata)
      !
      DO WHILE( xval>=Bkpt(ileft+1) .AND. ileft<n )
        ileft = ileft + 1
      END DO
      !
      CALL BSPLVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      irow = neqcon + np1 + nincon
      W(irow,1:n) = 0.E0
      intrvl = ileft - nordm1
      W(irow,intrvl:intrvl+Nord-1) = Bf(1:Nord,ideriv+1)
      !
      IF( itype==1 ) THEN
        W(irow,np1) = Yconst(idata)
      ELSE
        W(irow,np1) = -Yconst(idata)
        W(irow,intrvl:intrvl+Nord-1) = -W(irow,intrvl:intrvl+Nord-1)
      END IF
    END IF
  END DO
  !
  !     Solve constrained least squares equations.
  !
  CALL LSEI(W,Mdw,neqcon,np1,nincon,n,prgopt,Coeff,rnorme,rnorml,Mode,Work,Iwork)
END SUBROUTINE FCMN
