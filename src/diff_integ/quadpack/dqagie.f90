!** DQAGIE
PURE SUBROUTINE DQAGIE(F,Bound,Inf,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,&
    Alist,Blist,Rlist,Elist,Iord,Last)
  !> The routine calculates an approximation result to a given integral
  !  I = Integral of F over (BOUND,+INFINITY)
  !  or I = Integral of F over (-INFINITY,BOUND)
  !  or I = Integral of F over (-INFINITY,+INFINITY),
  !  hopefully satisfying following claim for accuracy
  !  ABS(I-RESULT)<=MAX(EPSABS,EPSREL*ABS(I))
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A3A1, H2A4A1
  !***
  ! **Type:**      DOUBLE PRECISION (QAGIE-S, DQAGIE-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
  !             GLOBALLY ADAPTIVE, INFINITE INTERVALS, QUADPACK,
  !             QUADRATURE, TRANSFORMATION
  !***
  ! **Author:**  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***
  ! **Description:**
  !
  ! Integration over infinite intervals
  ! Standard fortran subroutine
  !
  !            F      - Double precision
  !                     Function subprogram defining the integrand
  !                     function F(X). The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            BOUND  - Double precision
  !                     Finite bound of integration range
  !                     (has no meaning if interval is doubly-infinite)
  !
  !            INF    - Double precision
  !                     Indicating the kind of integration range involved
  !                     INF = 1 corresponds to  (BOUND,+INFINITY),
  !                     INF = -1            to  (-INFINITY,BOUND),
  !                     INF = 2             to (-INFINITY,+INFINITY).
  !
  !            EPSABS - Double precision
  !                     Absolute accuracy requested
  !            EPSREL - Double precision
  !                     Relative accuracy requested
  !                     If  EPSABS<=0
  !                     and EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     the routine will end with IER = 6.
  !
  !            LIMIT  - Integer
  !                     Gives an upper bound on the number of subintervals
  !                     in the partition of (A,B), LIMIT>=1
  !
  !         ON RETURN
  !            RESULT - Double precision
  !                     Approximation to the integral
  !
  !            ABSERR - Double precision
  !                     Estimate of the modulus of the absolute error,
  !                     which should equal or exceed ABS(I-RESULT)
  !
  !            NEVAL  - Integer
  !                     Number of integrand evaluations
  !
  !            IER    - Integer
  !                     IER = 0 Normal and reliable termination of the
  !                             routine. It is assumed that the requested
  !                             accuracy has been achieved.
  !                   - IER>0 Abnormal termination of the routine. The
  !                             estimates for result and error are less
  !                             reliable. It is assumed that the requested
  !                             accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more
  !                             subdivisions by increasing the value of
  !                             LIMIT (and taking the according dimension
  !                             adjustments into account).  However, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties.
  !                             If the position of a local difficulty can
  !                             be determined (e.g. SINGULARITY,
  !                             DISCONTINUITY within the interval) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling the
  !                             integrator on the subranges. If possible,
  !                             an appropriate special-purpose integrator
  !                             should be used, which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                             The error may be under-estimated.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 4 The algorithm does not converge.
  !                             Roundoff error is detected in the
  !                             extrapolation table.
  !                             It is assumed that the requested tolerance
  !                             cannot be achieved, and that the returned
  !                             result is the best which can be obtained.
  !                         = 5 The integral is probably divergent, or
  !                             slowly convergent. It must be noted that
  !                             divergence can occur with any other value
  !                             of IER.
  !                         = 6 The input is invalid, because
  !                             (EPSABS<=0 and
  !                              EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
  !                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
  !                             ELIST(1) and IORD(1) are set to zero.
  !                             ALIST(1) and BLIST(1) are set to 0
  !                             and 1 respectively.
  !
  !            ALIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the left
  !                     end points of the subintervals in the partition
  !                     of the transformed integration range (0,1).
  !
  !            BLIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the right
  !                     end points of the subintervals in the partition
  !                     of the transformed integration range (0,1).
  !
  !            RLIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the integral
  !                     approximations on the subintervals
  !
  !            ELIST  - Double precision
  !                     Vector of dimension at least LIMIT,  the first
  !                     LAST elements of which are the moduli of the
  !                     absolute error estimates on the subintervals
  !
  !            IORD   - Integer
  !                     Vector of dimension LIMIT, the first K
  !                     elements of which are pointers to the
  !                     error estimates over the subintervals,
  !                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
  !                     form a decreasing sequence, with K = LAST
  !                     If LAST<=(LIMIT/2+2), and K = LIMIT+1-LAST
  !                     otherwise
  !
  !            LAST   - Integer
  !                     Number of subintervals actually produced
  !                     in the subdivision process
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DQELG, DQK15I, DQPSRT

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  !
  INTERFACE
    REAL(DP) PURE FUNCTION F(X)
      IMPORT DP
      REAL(DP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Inf, Limit
  INTEGER, INTENT(OUT) :: Ier, Last, Neval, Iord(Limit)
  REAL(DP), INTENT(IN) :: Bound, Epsabs, Epsrel
  REAL(DP), INTENT(OUT) :: Abserr, Result
  REAL(DP), INTENT(OUT) :: Alist(Limit), Blist(Limit), Elist(Limit), Rlist(Limit)
  !
  INTEGER :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, maxerr, &
    nres, nrmax, numrl2
  REAL(DP) :: abseps, area, area1, area12, area2, a1, a2, boun, b1, b2, correc, &
    defabs, defab1, defab2, dres, epmach, erlarg, erlast, errbnd, errmax, error1, &
    error2, erro12, errsum, ertest, oflow, resabs, reseps, res3la(3), &
    rlist2(52), small, uflow
  LOGICAL :: extrap, noext
  !
  !            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
  !            LIMEXP IN SUBROUTINE DQELG.
  !
  !
  !            LIST OF MAJOR VARIABLES
  !            -----------------------
  !
  !           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
  !                       CONSIDERED UP TO NOW
  !           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
  !                       CONSIDERED UP TO NOW
  !           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
  !                       (ALIST(I),BLIST(I))
  !           RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+2),
  !                       CONTAINING THE PART OF THE EPSILON TABLE
  !                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
  !           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
  !           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
  !                       ESTIMATE
  !           ERRMAX    - ELIST(MAXERR)
  !           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
  !                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
  !           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
  !           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
  !           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
  !                       ABS(RESULT))
  !           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
  !           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
  !           LAST      - INDEX FOR SUBDIVISION
  !           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
  !           NUMRL2    - NUMBER OF ELEMENTS CURRENTLY IN RLIST2. IF AN
  !                       APPROPRIATE APPROXIMATION TO THE COMPOUNDED
  !                       INTEGRAL HAS BEEN OBTAINED, IT IS PUT IN
  !                       RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN INCREASED
  !                       BY ONE.
  !           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
  !                       TO NOW, MULTIPLIED BY 1.5
  !           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
  !                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
  !           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
  !                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
  !                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
  !                       TRY TO DECREASE THE VALUE OF ERLARG.
  !           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
  !                       IS NO LONGER ALLOWED (TRUE-VALUE)
  !
  !            MACHINE DEPENDENT CONSTANTS
  !            ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  DQAGIE
  epmach = D1MACH(4)
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           -----------------------------
  !
  Ier = 0
  Neval = 0
  Last = 0
  Result = 0._DP
  Abserr = 0._DP
  Alist(1) = 0._DP
  Blist(1) = 1._DP
  Rlist(1) = 0._DP
  Elist(1) = 0._DP
  Iord(1) = 0
  IF( Epsabs<=0._DP .AND. Epsrel<MAX(0.5E+02_DP*epmach,0.5E-28_DP) ) Ier = 6
  IF( Ier==6 ) RETURN
  !
  !
  !           FIRST APPROXIMATION TO THE INTEGRAL
  !           -----------------------------------
  !
  !           DETERMINE THE INTERVAL TO BE MAPPED ONTO (0,1).
  !           IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2, WHERE
  !           I1 = INTEGRAL OF F OVER (-INFINITY,0),
  !           I2 = INTEGRAL OF F OVER (0,+INFINITY).
  !
  boun = Bound
  IF( Inf==2 ) boun = 0._DP
  CALL DQK15I(F,boun,Inf,0._DP,1._DP,Result,Abserr,defabs,resabs)
  !
  !           TEST ON ACCURACY
  !
  Last = 1
  Rlist(1) = Result
  Elist(1) = Abserr
  Iord(1) = 1
  dres = ABS(Result)
  errbnd = MAX(Epsabs,Epsrel*dres)
  IF( Abserr<=100._DP*epmach*defabs .AND. Abserr>errbnd ) Ier = 2
  IF( Limit==1 ) Ier = 1
  IF( Ier/=0 .OR. (Abserr<=errbnd .AND. Abserr/=resabs) .OR. Abserr==0._DP )&
    GOTO 300
  !
  !           INITIALIZATION
  !           --------------
  !
  uflow = D1MACH(1)
  oflow = D1MACH(2)
  rlist2(1) = Result
  errmax = Abserr
  maxerr = 1
  area = Result
  errsum = Abserr
  Abserr = oflow
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .FALSE.
  noext = .FALSE.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  IF( dres>=(1._DP-0.5E+02_DP*epmach)*defabs ) ksgn = 1
  !
  !           MAIN DO-LOOP
  !           ------------
  !
  DO Last = 2, Limit
    !
    !           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
    !
    a1 = Alist(maxerr)
    b1 = 0.5_DP*(Alist(maxerr)+Blist(maxerr))
    a2 = b1
    b2 = Blist(maxerr)
    erlast = errmax
    CALL DQK15I(F,boun,Inf,a1,b1,area1,error1,resabs,defab1)
    CALL DQK15I(F,boun,Inf,a2,b2,area2,error2,resabs,defab2)
    !
    !           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
    !           AND ERROR AND TEST FOR ACCURACY.
    !
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - Rlist(maxerr)
    IF( defab1/=error1 .AND. defab2/=error2 ) THEN
      IF( ABS(Rlist(maxerr)-area12)<=0.1D-04*ABS(area12) .AND. &
          erro12>=0.99_DP*errmax ) THEN
        IF( extrap ) iroff2 = iroff2 + 1
        IF( .NOT. extrap ) iroff1 = iroff1 + 1
      END IF
      IF( Last>10 .AND. erro12>errmax ) iroff3 = iroff3 + 1
    END IF
    Rlist(maxerr) = area1
    Rlist(Last) = area2
    errbnd = MAX(Epsabs,Epsrel*ABS(area))
    !
    !           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
    !
    IF( iroff1+iroff2>=10 .OR. iroff3>=20 ) Ier = 2
    IF( iroff2>=5 ) ierro = 3
    !
    !           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
    !           SUBINTERVALS EQUALS LIMIT.
    !
    IF( Last==Limit ) Ier = 1
    !
    !           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
    !           AT SOME POINTS OF THE INTEGRATION RANGE.
    !
    IF( MAX(ABS(a1),ABS(b2))<=(1._DP+100._DP*epmach)&
      *(ABS(a2)+0.1D+04*uflow) ) Ier = 4
    !
    !           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
    !
    IF( error2>error1 ) THEN
      Alist(maxerr) = a2
      Alist(Last) = a1
      Blist(Last) = b1
      Rlist(maxerr) = area2
      Rlist(Last) = area1
      Elist(maxerr) = error2
      Elist(Last) = error1
    ELSE
      Alist(Last) = a2
      Blist(maxerr) = b1
      Blist(Last) = b2
      Elist(maxerr) = error1
      Elist(Last) = error2
    END IF
    !
    !           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
    !           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
    !           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
    !
    CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
    IF( errsum<=errbnd ) GOTO 200
    IF( Ier/=0 ) EXIT
    IF( Last==2 ) THEN
      small = 0.375_DP
      erlarg = errsum
      ertest = errbnd
      rlist2(2) = area
    ELSEIF( .NOT. (noext) ) THEN
      erlarg = erlarg - erlast
      IF( ABS(b1-a1)>small ) erlarg = erlarg + erro12
      IF( .NOT. (extrap) ) THEN
        !
        !           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
        !           SMALLEST INTERVAL.
        !
        IF( ABS(Blist(maxerr)-Alist(maxerr))>small ) CYCLE
        extrap = .TRUE.
        nrmax = 2
      END IF
      IF( ierro/=3 .AND. erlarg>ertest ) THEN
        !
        !           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
        !           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
        !           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
        !
        id = nrmax
        jupbnd = Last
        IF( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
        DO k = id, jupbnd
          maxerr = Iord(nrmax)
          errmax = Elist(maxerr)
          IF( ABS(Blist(maxerr)-Alist(maxerr))>small ) GOTO 100
          nrmax = nrmax + 1
        END DO
      END IF
      !
      !           PERFORM EXTRAPOLATION.
      !
      numrl2 = numrl2 + 1
      rlist2(numrl2) = area
      CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
      ktmin = ktmin + 1
      IF( ktmin>5 .AND. Abserr<0.1D-02*errsum ) Ier = 5
      IF( abseps<Abserr ) THEN
        ktmin = 0
        Abserr = abseps
        Result = reseps
        correc = erlarg
        ertest = MAX(Epsabs,Epsrel*ABS(reseps))
        IF( Abserr<=ertest ) EXIT
      END IF
      !
      !            PREPARE BISECTION OF THE SMALLEST INTERVAL.
      !
      IF( numrl2==1 ) noext = .TRUE.
      IF( Ier==5 ) EXIT
      maxerr = Iord(1)
      errmax = Elist(maxerr)
      nrmax = 1
      extrap = .FALSE.
      small = small*0.5_DP
      erlarg = errsum
    END IF
    100 CONTINUE
  END DO
  !
  !           SET FINAL RESULT AND ERROR ESTIMATE.
  !           ------------------------------------
  !
  IF( Abserr/=oflow ) THEN
    IF( (Ier+ierro)/=0 ) THEN
      IF( ierro==3 ) Abserr = Abserr + correc
      IF( Ier==0 ) Ier = 3
      IF( Result==0._DP .OR. area==0._DP ) THEN
        IF( Abserr>errsum ) GOTO 200
        IF( area==0._DP ) GOTO 300
      ELSEIF( Abserr/ABS(Result)>errsum/ABS(area) ) THEN
        GOTO 200
      END IF
    END IF
    !
    !           TEST ON DIVERGENCE
    !
    IF( ksgn/=(-1) .OR. MAX(ABS(Result),ABS(area))>defabs*0.1D-01 ) THEN
      IF( 0.1D-01>(Result/area) .OR. (Result/area)>0.1D+03 .OR. &
        errsum>ABS(area) ) Ier = 6
    END IF
    GOTO 300
  END IF
  !
  !           COMPUTE GLOBAL INTEGRAL SUM.
  !
  200  Result = 0._DP
  DO k = 1, Last
    Result = Result + Rlist(k)
  END DO
  Abserr = errsum
  300  Neval = 30*Last - 15
  IF( Inf==2 ) Neval = 2*Neval
  IF( Ier>2 ) Ier = Ier - 1
  !
  RETURN
END SUBROUTINE DQAGIE