!DECK DDSTP
SUBROUTINE DDSTP(Eps,F,FA,Hmax,Impl,Ierror,JACOBN,Matdim,Maxord,Mint,&
    Miter,Ml,Mu,N,Nde,Ywt,Uround,USERS,Avgh,Avgord,H,Hused,&
    Jtask,Mntold,Mtrold,Nfe,Nje,Nqused,Nstep,T,Y,Yh,A,Convrg,&
    Dfdy,El,Fac,Hold,Ipvt,Jstate,Jstepl,Nq,Nwait,Rc,Rmax,&
    Save1,Save2,Tq,Trend,Iswflg,Mtrsv,Mxrdsv)
  IMPLICIT NONE
  REAL FA, USERS
  INTEGER JACOBN
  !***BEGIN PROLOGUE  DDSTP
  !***SUBSIDIARY
  !***PURPOSE  DDSTP performs one step of the integration of an initial
  !            value problem for a system of ordinary differential
  !            equations.
  !***LIBRARY   SLATEC (SDRIVE)
  !***TYPE      DOUBLE PRECISION (SDSTP-S, DDSTP-D, CDSTP-C)
  !***AUTHOR  Kahaner, D. K., (NIST)
  !             National Institute of Standards and Technology
  !             Gaithersburg, MD  20899
  !           Sutherland, C. D., (LANL)
  !             Mail Stop D466
  !             Los Alamos National Laboratory
  !             Los Alamos, NM  87545
  !***DESCRIPTION
  !
  !  Communication with DDSTP is done with the following variables:
  !
  !    YH      An N by MAXORD+1 array containing the dependent variables
  !              and their scaled derivatives.  MAXORD, the maximum order
  !              used, is currently 12 for the Adams methods and 5 for the
  !              Gear methods.  YH(I,J+1) contains the J-th derivative of
  !              Y(I), scaled by H**J/factorial(J).  Only Y(I),
  !              1 .LE. I .LE. N, need be set by the calling program on
  !              the first entry.  The YH array should not be altered by
  !              the calling program.  When referencing YH as a
  !              2-dimensional array, use a column length of N, as this is
  !              the value used in DDSTP.
  !    DFDY    A block of locations used for partial derivatives if MITER
  !              is not 0.  If MITER is 1 or 2 its length must be at least
  !              N*N.  If MITER is 4 or 5 its length must be at least
  !              (2*ML+MU+1)*N.
  !    YWT     An array of N locations used in convergence and error tests
  !    SAVE1
  !    SAVE2   Arrays of length N used for temporary storage.
  !    IPVT    An integer array of length N used by the linear system
  !              solvers for the storage of row interchange information.
  !    A       A block of locations used to store the matrix A, when using
  !              the implicit method.  If IMPL is 1, A is a MATDIM by N
  !              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4
  !              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N.
  !              If IMPL is 3, A is a MATDIM by NDE array.
  !    JTASK   An integer used on input.
  !              It has the following values and meanings:
  !                 .EQ. 0  Perform the first step.  This value enables
  !                         the subroutine to initialize itself.
  !                .GT. 0  Take a new step continuing from the last.
  !                         Assumes the last step was successful and
  !                         user has not changed any parameters.
  !                 .LT. 0  Take a new step with a new value of H and/or
  !                         MINT and/or MITER.
  !    JSTATE  A completion code with the following meanings:
  !                1  The step was successful.
  !                2  A solution could not be obtained with H .NE. 0.
  !                3  A solution was not obtained in MXTRY attempts.
  !                4  For IMPL .NE. 0, the matrix A is singular.
  !              On a return with JSTATE .GT. 1, the values of T and
  !              the YH array are as of the beginning of the last
  !              step, and H is the last step size attempted.
  !
  !***ROUTINES CALLED  DDCOR, DDCST, DDNTL, DDPSC, DDPST, DDSCL, DNRM2
  !***REVISION HISTORY  (YYMMDD)
  !   790601  DATE WRITTEN
  !   900329  Initial submission to SLATEC.
  !***END PROLOGUE  DDSTP
  EXTERNAL F, JACOBN, FA, USERS
  INTEGER i, Ierror, Impl, Ipvt(*), Iswflg, iter, j, Jstate, &
    Jstepl, Jtask, Matdim, Maxord, Mint, Miter, Ml, Mntold, &
    Mtrold, Mtrsv, Mu, MXFAIL, MXITER, Mxrdsv, MXTRY, N, Nde, &
    NDJSTP, nfail, Nfe, Nje, Nq, Nqused, Nstep, nsv, ntry, &
    Nwait
  REAL(8) :: A(Matdim,*), Avgh, Avgord, BIAS1, BIAS2, BIAS3, &
    bnd, ctest, d, denom, Dfdy(Matdim,*), d1, El(13,12)&
    , Eps, erdn, erup, etest, Fac(*), H, Hmax, hn, &
    Hold, hs, Hused, numer, Rc, RCTEST, rh, rh1, &
    rh2, rh3, Rmax, RMFAIL, RMNORM, Save1(*), Save2(*)&
    , DNRM2, T, told, Tq(3,12), Trend, TRSHLD, Uround, &
    Y(*), Yh(N,*), Ywt(*), y0nrm
  LOGICAL Convrg, evalfa, evaljc, ier, switch
  PARAMETER (BIAS1=1.3D0,BIAS2=1.2D0,BIAS3=1.4D0,MXFAIL=3,MXITER=3,MXTRY=50,&
    RCTEST=.3D0,RMFAIL=2.D0,RMNORM=10.D0,TRSHLD=1.D0)
  PARAMETER (NDJSTP=10)
  DATA ier/.FALSE./
  !***FIRST EXECUTABLE STATEMENT  DDSTP
  nsv = N
  bnd = 0.D0
  switch = .FALSE.
  ntry = 0
  told = T
  nfail = 0
  IF ( Jtask<=0 ) THEN
    CALL DDNTL(Eps,F,FA,Hmax,Hold,Impl,Jtask,Matdim,Maxord,Mint,Miter,Ml,Mu,&
      N,Nde,Save1,T,Uround,USERS,Y,Ywt,H,Mntold,Mtrold,Nfe,Rc,Yh,A,&
      Convrg,El,Fac,ier,Ipvt,Nq,Nwait,rh,Rmax,Save2,Tq,Trend,&
      Iswflg,Jstate)
    IF ( N==0 ) GOTO 800
    IF ( H==0.D0 ) GOTO 500
    IF ( ier ) GOTO 600
  ENDIF
  100  ntry = ntry + 1
  IF ( ntry>MXTRY ) THEN
    !
    Jstate = 3
    Hold = H
    RETURN
  ELSE
    T = T + H
    CALL DDPSC(1,N,Nq,Yh)
    evaljc = (((ABS(Rc-1.D0)>RCTEST).OR.(Nstep>=Jstepl+NDJSTP)).AND.&
      (Miter/=0))
    evalfa = .NOT.evaljc
  ENDIF
  !
  200  iter = 0
  DO i = 1, N
    Y(i) = Yh(i,1)
  ENDDO
  CALL F(N,T,Y,Save2)
  IF ( N==0 ) THEN
    Jstate = 6
    GOTO 700
  ENDIF
  Nfe = Nfe + 1
  IF ( evaljc.OR.ier ) THEN
    CALL DDPST(El,F,FA,H,Impl,JACOBN,Matdim,Miter,Ml,Mu,N,Nde,Nq,Save2,T,&
      USERS,Y,Yh,Ywt,Uround,Nfe,Nje,A,Dfdy,Fac,ier,Ipvt,Save1,&
      Iswflg,bnd,Jstate)
    IF ( N==0 ) GOTO 700
    IF ( ier ) GOTO 300
    Convrg = .FALSE.
    Rc = 1.D0
    Jstepl = Nstep
  ENDIF
  DO i = 1, N
    Save1(i) = 0.D0
  ENDDO
  DO
    !                      Up to MXITER corrector iterations are taken.
    !                      Convergence is tested by requiring the r.m.s.
    !                      norm of changes to be less than EPS.  The sum of
    !                      the corrections is accumulated in the vector
    !                      SAVE1(I).  It is approximately equal to the L-th
    !                      derivative of Y multiplied by
    !                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
    !                      proportional to the actual errors to the lowest
    !                      power of H present (H**L).  The YH array is not
    !                      altered in the correction loop.  The norm of the
    !                      iterate difference is stored in D.  If
    !                      ITER .GT. 0, an estimate of the convergence rate
    !                      constant is stored in TREND, and this is used in
    !                      the convergence test.
    !
    CALL DDCOR(Dfdy,El,FA,H,Ierror,Impl,Ipvt,Matdim,Miter,Ml,Mu,N,Nde,Nq,T,&
      USERS,Y,Yh,Ywt,evalfa,Save1,Save2,A,d,Jstate)
    IF ( N==0 ) GOTO 700
    IF ( Iswflg==3.AND.Mint==1 ) THEN
      IF ( iter==0 ) THEN
        numer = DNRM2(N,Save1,1)
        DO i = 1, N
          Dfdy(1,i) = Save1(i)
        ENDDO
        y0nrm = DNRM2(N,Yh,1)
      ELSE
        denom = numer
        DO i = 1, N
          Dfdy(1,i) = Save1(i) - Dfdy(1,i)
        ENDDO
        numer = DNRM2(N,Dfdy,Matdim)
        IF ( El(1,Nq)*numer<=100.D0*Uround*y0nrm ) THEN
          IF ( Rmax==RMFAIL ) THEN
            switch = .TRUE.
            GOTO 400
          ENDIF
        ENDIF
        DO i = 1, N
          Dfdy(1,i) = Save1(i)
        ENDDO
        IF ( denom/=0.D0 ) bnd = MAX(bnd,numer/(denom*ABS(H)*El(1,Nq)))
      ENDIF
    ENDIF
    IF ( iter>0 ) Trend = MAX(.9D0*Trend,d/d1)
    d1 = d
    ctest = MIN(2.D0*Trend,1.D0)*d
    IF ( ctest<=Eps ) GOTO 400
    iter = iter + 1
    IF ( iter<MXITER ) THEN
      DO i = 1, N
        Y(i) = Yh(i,1) + El(1,Nq)*Save1(i)
      ENDDO
      CALL F(N,T,Y,Save2)
      IF ( N==0 ) THEN
        Jstate = 6
        GOTO 700
      ENDIF
      Nfe = Nfe + 1
      CYCLE
    ENDIF
    !                     The corrector iteration failed to converge in
    !                     MXITER tries.  If partials are involved but are
    !                     not up to date, they are reevaluated for the next
    !                     try.  Otherwise the YH array is retracted to its
    !                     values before prediction, and H is reduced, if
    !                     possible.  If not, a no-convergence exit is taken.
    IF ( Convrg ) THEN
      evaljc = .TRUE.
      evalfa = .FALSE.
      GOTO 200
    ENDIF
    EXIT
  ENDDO
  300  T = told
  CALL DDPSC(-1,N,Nq,Yh)
  Nwait = Nq + 2
  IF ( Jtask/=0.AND.Jtask/=2 ) Rmax = RMFAIL
  IF ( iter==0 ) THEN
    rh = .3D0
  ELSE
    rh = .9D0*(Eps/ctest)**(.2D0)
  ENDIF
  IF ( rh*H==0.D0 ) GOTO 500
  CALL DDSCL(Hmax,N,Nq,Rmax,H,Rc,rh,Yh)
  GOTO 100
  !                          The corrector has converged.  CONVRG is set
  !                          to .TRUE. if partial derivatives were used,
  !                          to indicate that they may need updating on
  !                          subsequent steps.  The error test is made.
  400  Convrg = (Miter/=0)
  IF ( Ierror==1.OR.Ierror==5 ) THEN
    DO i = 1, Nde
      Save2(i) = Save1(i)/Ywt(i)
    ENDDO
  ELSE
    DO i = 1, Nde
      Save2(i) = Save1(i)/MAX(ABS(Y(i)),Ywt(i))
    ENDDO
  ENDIF
  etest = DNRM2(Nde,Save2,1)/(Tq(2,Nq)*SQRT(REAL(Nde, 8)))
  !
  !                           The error test failed.  NFAIL keeps track of
  !                           multiple failures.  Restore T and the YH
  !                           array to their previous values, and prepare
  !                           to try the step again.  Compute the optimum
  !                           step size for this or one lower order.
  IF ( etest>Eps ) THEN
    T = told
    CALL DDPSC(-1,N,Nq,Yh)
    nfail = nfail + 1
    IF ( nfail<MXFAIL.OR.Nq==1 ) THEN
      IF ( Jtask/=0.AND.Jtask/=2 ) Rmax = RMFAIL
      rh2 = 1.D0/(BIAS2*(etest/Eps)**(1.D0/(Nq+1)))
      IF ( Nq>1 ) THEN
        IF ( Ierror==1.OR.Ierror==5 ) THEN
          DO i = 1, Nde
            Save2(i) = Yh(i,Nq+1)/Ywt(i)
          ENDDO
        ELSE
          DO i = 1, Nde
            Save2(i) = Yh(i,Nq+1)/MAX(ABS(Y(i)),Ywt(i))
          ENDDO
        ENDIF
        erdn = DNRM2(Nde,Save2,1)/(Tq(1,Nq)*SQRT(REAL(Nde, 8)))
        rh1 = 1.D0/MAX(1.D0,BIAS1*(erdn/Eps)**(1.D0/Nq))
        IF ( rh2<rh1 ) THEN
          Nq = Nq - 1
          Rc = Rc*El(1,Nq)/El(1,Nq+1)
          rh = rh1
        ELSE
          rh = rh2
        ENDIF
      ELSE
        rh = rh2
      ENDIF
      Nwait = Nq + 2
      IF ( rh*H==0.D0 ) GOTO 500
      CALL DDSCL(Hmax,N,Nq,Rmax,H,Rc,rh,Yh)
      GOTO 100
    ENDIF
    !                Control reaches this section if the error test has
    !                failed MXFAIL or more times.  It is assumed that the
    !                derivatives that have accumulated in the YH array have
    !                errors of the wrong order.  Hence the first derivative
    !                is recomputed, the order is set to 1, and the step is
    !                retried.
    nfail = 0
    Jtask = 2
    DO i = 1, N
      Y(i) = Yh(i,1)
    ENDDO
    CALL DDNTL(Eps,F,FA,Hmax,Hold,Impl,Jtask,Matdim,Maxord,Mint,Miter,Ml,Mu,&
      N,Nde,Save1,T,Uround,USERS,Y,Ywt,H,Mntold,Mtrold,Nfe,Rc,Yh,A,&
      Convrg,El,Fac,ier,Ipvt,Nq,Nwait,rh,Rmax,Save2,Tq,Trend,&
      Iswflg,Jstate)
    Rmax = RMNORM
    IF ( N==0 ) GOTO 800
    IF ( H==0.D0 ) GOTO 500
    IF ( .NOT.(ier) ) GOTO 100
    GOTO 600
  ENDIF
  !                          After a successful step, update the YH array.
  Nstep = Nstep + 1
  Hused = H
  Nqused = Nq
  Avgh = ((Nstep-1)*Avgh+H)/Nstep
  Avgord = ((Nstep-1)*Avgord+Nq)/Nstep
  DO j = 1, Nq + 1
    DO i = 1, N
      Yh(i,j) = Yh(i,j) + El(j,Nq)*Save1(i)
    ENDDO
  ENDDO
  DO i = 1, N
    Y(i) = Yh(i,1)
  ENDDO
  !                                          If ISWFLG is 3, consider
  !                                          changing integration methods.
  IF ( Iswflg==3 ) THEN
    IF ( bnd/=0.D0 ) THEN
      IF ( Mint==1.AND.Nq<=5 ) THEN
        hn = ABS(H)/MAX(Uround,(etest/Eps)**(1.D0/(Nq+1)))
        hn = MIN(hn,1.D0/(2.D0*El(1,Nq)*bnd))
        hs = ABS(H)/MAX(Uround,(etest/(Eps*El(Nq+1,1)))**(1.D0/(Nq+1)))
        IF ( hs>1.2D0*hn ) THEN
          Mint = 2
          Mntold = Mint
          Miter = Mtrsv
          Mtrold = Miter
          Maxord = MIN(Mxrdsv,5)
          Rc = 0.D0
          Rmax = RMNORM
          Trend = 1.D0
          CALL DDCST(Maxord,Mint,Iswflg,El,Tq)
          Nwait = Nq + 2
        ENDIF
      ELSEIF ( Mint==2 ) THEN
        hs = ABS(H)/MAX(Uround,(etest/Eps)**(1.D0/(Nq+1)))
        hn = ABS(H)/MAX(Uround,(etest*El(Nq+1,1)/Eps)**(1.D0/(Nq+1)))
        hn = MIN(hn,1.D0/(2.D0*El(1,Nq)*bnd))
        IF ( hn>=hs ) THEN
          Mint = 1
          Mntold = Mint
          Miter = 0
          Mtrold = Miter
          Maxord = MIN(Mxrdsv,12)
          Rmax = RMNORM
          Trend = 1.D0
          Convrg = .FALSE.
          CALL DDCST(Maxord,Mint,Iswflg,El,Tq)
          Nwait = Nq + 2
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  IF ( switch ) THEN
    Mint = 2
    Mntold = Mint
    Miter = Mtrsv
    Mtrold = Miter
    Maxord = MIN(Mxrdsv,5)
    Nq = MIN(Nq,Maxord)
    Rc = 0.D0
    Rmax = RMNORM
    Trend = 1.D0
    CALL DDCST(Maxord,Mint,Iswflg,El,Tq)
    Nwait = Nq + 2
  ENDIF
  !                           Consider changing H if NWAIT = 1.  Otherwise
  !                           decrease NWAIT by 1.  If NWAIT is then 1 and
  !                           NQ.LT.MAXORD, then SAVE1 is saved for use in
  !                           a possible order increase on the next step.
  !
  IF ( Jtask==0.OR.Jtask==2 ) THEN
    rh = 1.D0/MAX(Uround,BIAS2*(etest/Eps)**(1.D0/(Nq+1)))
    IF ( rh>TRSHLD ) CALL DDSCL(Hmax,N,Nq,Rmax,H,Rc,rh,Yh)
  ELSEIF ( Nwait>1 ) THEN
    Nwait = Nwait - 1
    IF ( Nwait==1.AND.Nq<Maxord ) THEN
      DO i = 1, Nde
        Yh(i,Maxord+1) = Save1(i)
      ENDDO
    ENDIF
    !             If a change in H is considered, an increase or decrease in
    !             order by one is considered also.  A change in H is made
    !             only if it is by a factor of at least TRSHLD.  Factors
    !             RH1, RH2, and RH3 are computed, by which H could be
    !             multiplied at order NQ - 1, order NQ, or order NQ + 1,
    !             respectively.  The largest of these is determined and the
    !             new order chosen accordingly.  If the order is to be
    !             increased, we compute one additional scaled derivative.
    !             If there is a change of order, reset NQ and the
    !             coefficients.  In any case H is reset according to RH and
    !             the YH array is rescaled.
  ELSE
    IF ( Nq==1 ) THEN
      rh1 = 0.D0
    ELSE
      IF ( Ierror==1.OR.Ierror==5 ) THEN
        DO i = 1, Nde
          Save2(i) = Yh(i,Nq+1)/Ywt(i)
        ENDDO
      ELSE
        DO i = 1, Nde
          Save2(i) = Yh(i,Nq+1)/MAX(ABS(Y(i)),Ywt(i))
        ENDDO
      ENDIF
      erdn = DNRM2(Nde,Save2,1)/(Tq(1,Nq)*SQRT(REAL(Nde, 8)))
      rh1 = 1.D0/MAX(Uround,BIAS1*(erdn/Eps)**(1.D0/Nq))
    ENDIF
    rh2 = 1.D0/MAX(Uround,BIAS2*(etest/Eps)**(1.D0/(Nq+1)))
    IF ( Nq==Maxord ) THEN
      rh3 = 0.D0
    ELSE
      IF ( Ierror==1.OR.Ierror==5 ) THEN
        DO i = 1, Nde
          Save2(i) = (Save1(i)-Yh(i,Maxord+1))/Ywt(i)
        ENDDO
      ELSE
        DO i = 1, Nde
          Save2(i) = (Save1(i)-Yh(i,Maxord+1))/MAX(ABS(Y(i)),Ywt(i))
        ENDDO
      ENDIF
      erup = DNRM2(Nde,Save2,1)/(Tq(3,Nq)*SQRT(REAL(Nde, 8)))
      rh3 = 1.D0/MAX(Uround,BIAS3*(erup/Eps)**(1.D0/(Nq+2)))
    ENDIF
    IF ( rh1>rh2.AND.rh1>=rh3 ) THEN
      rh = rh1
      IF ( rh<=TRSHLD ) GOTO 450
      Nq = Nq - 1
      Rc = Rc*El(1,Nq)/El(1,Nq+1)
    ELSEIF ( rh2>=rh1.AND.rh2>=rh3 ) THEN
      rh = rh2
      IF ( rh<=TRSHLD ) GOTO 450
    ELSE
      rh = rh3
      IF ( rh<=TRSHLD ) GOTO 450
      DO i = 1, N
        Yh(i,Nq+2) = Save1(i)*El(Nq+1,Nq)/(Nq+1)
      ENDDO
      Nq = Nq + 1
      Rc = Rc*El(1,Nq)/El(1,Nq-1)
    ENDIF
    IF ( Iswflg==3.AND.Mint==1 ) THEN
      IF ( bnd/=0.D0 ) rh = MIN(rh,1.D0/(2.D0*El(1,Nq)*bnd*ABS(H)))
    ENDIF
    CALL DDSCL(Hmax,N,Nq,Rmax,H,Rc,rh,Yh)
    Rmax = RMNORM
    450    Nwait = Nq + 2
  ENDIF
  !               All returns are made through this section.  H is saved
  !               in HOLD to allow the caller to change H on the next step
  Jstate = 1
  Hold = H
  RETURN
  !
  500  Jstate = 2
  Hold = H
  DO i = 1, N
    Y(i) = Yh(i,1)
  ENDDO
  RETURN
  !
  600  Jstate = 4
  Hold = H
  RETURN
  !
  700  T = told
  CALL DDPSC(-1,nsv,Nq,Yh)
  DO i = 1, nsv
    Y(i) = Yh(i,1)
  ENDDO
  800  Hold = H
END SUBROUTINE DDSTP
