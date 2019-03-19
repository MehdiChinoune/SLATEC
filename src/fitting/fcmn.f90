!** FCMN
SUBROUTINE FCMN(Ndata,Xdata,Ydata,Sddata,Nord,Nbkpt,Bkptin,Nconst,Xconst,&
    Yconst,Nderiv,Mode,Coeff,Bf,Xtemp,Ptemp,Bkpt,G,Mdg,W,Mdw,&
    Work,Iwork)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to FC
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
  
  INTEGER Iwork(*), Mdg, Mdw, Mode, Nbkpt, Nconst, Ndata, Nderiv(*), &
    Nord
  REAL Bf(Nord,*), Bkpt(*), Bkptin(*), Coeff(*), G(Mdg,*), Ptemp(*), &
    Sddata(*), W(Mdw,*), Work(*), Xconst(*), Xdata(*), Xtemp(*), &
    Yconst(*), Ydata(*)
  !
  EXTERNAL BNDACC, BNDSOL, BSPLVD, BSPLVN, LSEI, SAXPY, SCOPY, &
    SSCAL, SSORT, XERMSG
  !
  REAL dummy, prgopt(10), rnorm, rnorme, rnorml, xmax, xmin, xval, &
    yval
  INTEGER i, idata, ideriv, ileft, intrvl, intw1, ip, ir, irow, &
    itype, iw1, iw2, l, lw, mt, n, nb, neqcon, nincon, &
    nordm1, nordp1, np1
  LOGICAL band, new, var
  CHARACTER(8) :: xern1
  !
  !* FIRST EXECUTABLE STATEMENT  FCMN
  !
  !     Analyze input.
  !
  IF ( Nord<1.OR.Nord>20 ) THEN
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, THE ORDER OF THE B-SPLINE MUST BE 1 THRU 20.',2,1)
    Mode = -1
    RETURN
    !
  ELSEIF ( Nbkpt<2*Nord ) THEN
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, THE NUMBER OF KNOTS MUST BE AT LEAST TWICE THE B-SPLINE ORDER.',2,1)
    Mode = -1
    RETURN
  ENDIF
  !
  IF ( Ndata<0 ) THEN
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, THE NUMBER OF DATA POINTS MUST BE NONNEGATIVE.',2,1)
    Mode = -1
    RETURN
  ENDIF
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
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK NB = '//xern1,&
      2,1)
    Mode = -1
    RETURN
  ENDIF
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
    CALL XERMSG('SLATEC','FCMN','IN FC, INPUT VALUE OF MODE MUST BE 1-4.',2,&
      1)
    Mode = -1
    RETURN
  ENDIF
  Mode = 0
  !
  !     Sort the breakpoints.
  !
  CALL SCOPY(Nbkpt,Bkptin,1,Bkpt,1)
  CALL SSORT(Bkpt,[dummy],Nbkpt,1)
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
    ENDIF
  ENDDO
  !
  !     Compute the number of variables.
  !
  n = Nbkpt - Nord
  np1 = n + 1
  lw = nb + (np1+Nconst)*np1 + 2*(neqcon+np1) + (nincon+np1) + (nincon+2)&
    *(np1+6)
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
  ENDDO
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
  IF ( var ) THEN
    prgopt(3) = 1
  ELSE
    prgopt(3) = 0
  ENDIF
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
  IF ( new ) THEN
    !
    !        To process least squares equations sort data and an array of
    !        pointers.
    !
    CALL SCOPY(Ndata,Xdata,1,Xtemp,1)
    DO i = 1, Ndata
      Ptemp(i) = i
    ENDDO
    !
    IF ( Ndata>0 ) THEN
      CALL SSORT(Xtemp,Ptemp,Ndata,2)
      xmin = MIN(xmin,Xtemp(1))
      xmax = MAX(xmax,Xtemp(Ndata))
    ENDIF
    !
    !        Fix breakpoint array if needed.
    !
    DO i = 1, Nord
      Bkpt(i) = MIN(Bkpt(i),xmin)
    ENDDO
    !
    DO i = np1, Nbkpt
      Bkpt(i) = MAX(Bkpt(i),xmax)
    ENDDO
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
      IF ( xval>=Bkpt(ileft+1) ) THEN
        CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
        !
        !              Move pointer up to have BKPT(ILEFT).LE.XVAL,
        !                 ILEFT.LT.NP1.
        !
        DO WHILE ( xval>=Bkpt(ileft+1).AND.ileft<n )
          ileft = ileft + 1
        ENDDO
      ENDIF
      !
      !           Obtain B-spline function value.
      !
      CALL BSPLVN(Bkpt,Nord,1,xval,ileft,Bf)
      !
      !           Move row into place.
      !
      irow = ir + mt
      mt = mt + 1
      CALL SCOPY(Nord,Bf,1,G(irow,1),Mdg)
      G(irow,nordp1) = Ydata(l)
      !
      !           Scale data if uncertainty is nonzero.
      !
      IF ( Sddata(l)/=0.E0 ) CALL SSCAL(nordp1,1.E0/Sddata(l),G(irow,1),Mdg)
      !
      !           When staging work area is exhausted, process rows.
      !
      IF ( irow==Mdg-1 ) THEN
        CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
        mt = 0
      ENDIF
    ENDDO
    !
    !        Process last block of equations.
    !
    CALL BNDACC(G,Mdg,Nord,ip,ir,mt,ileft-nordm1)
    !
    !        Last call to adjust block positioning.
    !
    CALL SCOPY(nordp1,0.E0,0,G(ir,1),Mdg)
    CALL BNDACC(G,Mdg,Nord,ip,ir,1,np1)
  ENDIF
  !
  band = band .AND. Nconst==0
  DO i = 1, n
    band = band .AND. G(i,1)/=0.E0
  ENDDO
  !
  !     Process banded least squares equations.
  !
  IF ( band ) THEN
    CALL BNDSOL(1,G,Mdg,Nord,ip,ir,Coeff,n,rnorm)
    RETURN
  ENDIF
  !
  !     Check further for sufficient storage in working arrays.
  !
  IF ( iw1<lw ) THEN
    WRITE (xern1,'(I8)') lw
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR W(*).  CHECK LW = '//xern1,&
      2,1)
    Mode = -1
    RETURN
  ENDIF
  !
  IF ( iw2<intw1 ) THEN
    WRITE (xern1,'(I8)') intw1
    CALL XERMSG('SLATEC','FCMN',&
      'IN FC, INSUFFICIENT STORAGE FOR IW(*).  CHECK IW1 = '//&
      xern1,2,1)
    Mode = -1
    RETURN
  ENDIF
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
      ENDDO
      !
      CALL BSPLVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      CALL SCOPY(np1,0.E0,0,W(neqcon,1),Mdw)
      CALL SCOPY(Nord,Bf(1,ideriv+1),1,W(neqcon,ileft-nordm1),Mdw)
      !
      IF ( itype==2 ) THEN
        W(neqcon,np1) = Yconst(idata)
      ELSE
        ileft = Nord
        yval = Yconst(idata)
        !
        DO WHILE ( yval>=Bkpt(ileft+1).AND.ileft<n )
          ileft = ileft + 1
        ENDDO
        !
        CALL BSPLVD(Bkpt,Nord,yval,ileft,Bf,ideriv+1)
        CALL SAXPY(Nord,-1.E0,Bf(1,ideriv+1),1,W(neqcon,ileft-nordm1),Mdw)
      ENDIF
    ENDIF
  ENDDO
  !
  !     Transfer least squares data.
  !
  DO i = 1, np1
    irow = i + neqcon
    CALL SCOPY(n,0.E0,0,W(irow,1),Mdw)
    CALL SCOPY(MIN(np1-i,Nord),G(i,1),Mdg,W(irow,i),Mdw)
    W(irow,np1) = G(i,nordp1)
  ENDDO
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
      ENDDO
      !
      CALL BSPLVD(Bkpt,Nord,xval,ileft,Bf,ideriv+1)
      irow = neqcon + np1 + nincon
      CALL SCOPY(n,0.E0,0,W(irow,1),Mdw)
      intrvl = ileft - nordm1
      CALL SCOPY(Nord,Bf(1,ideriv+1),1,W(irow,intrvl),Mdw)
      !
      IF ( itype==1 ) THEN
        W(irow,np1) = Yconst(idata)
      ELSE
        W(irow,np1) = -Yconst(idata)
        CALL SSCAL(Nord,-1.E0,W(irow,intrvl),Mdw)
      ENDIF
    ENDIF
  ENDDO
  !
  !     Solve constrained least squares equations.
  !
  CALL LSEI(W,Mdw,neqcon,np1,nincon,n,prgopt,Coeff,rnorme,rnorml,Mode,Work,&
    Iwork)
END SUBROUTINE FCMN
