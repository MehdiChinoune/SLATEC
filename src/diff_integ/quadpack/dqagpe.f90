!** DQAGPE
SUBROUTINE DQAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,Limit,Result,Abserr,&
    Neval,Ier,Alist,Blist,Rlist,Elist,Pts,Iord,Level,Ndin,Last)
  !> Approximate a given definite integral I = Integral of F
  !            over (A,B), hopefully satisfying the accuracy claim:
  !                 ABS(I-RESULT)<=MAX(EPSABS,EPSREL*ABS(I)).
  !            Break points of the integration interval, where local
  !            difficulties of the integrand may occur (e.g. singularities
  !            or discontinuities) are provided by the user.
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A1
  !***
  ! **Type:**      DOUBLE PRECISION (QAGPE-S, DQAGPE-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
  !             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
  !             SINGULARITIES AT USER SPECIFIED POINTS
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
  !        Computation of a definite integral
  !        Standard fortran subroutine
  !        Double precision version
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Double precision
  !                     Function subprogram defining the integrand
  !                     function F(X). The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            A      - Double precision
  !                     Lower limit of integration
  !
  !            B      - Double precision
  !                     Upper limit of integration
  !
  !            NPTS2  - Integer
  !                     Number equal to two more than the number of
  !                     user-supplied break points within the integration
  !                     range, NPTS2>=2.
  !                     If NPTS2<2, the routine will end with IER = 6.
  !
  !            POINTS - Double precision
  !                     Vector of dimension NPTS2, the first (NPTS2-2)
  !                     elements of which are the user provided break
  !                     POINTS. If these POINTS do not constitute an
  !                     ascending sequence there will be an automatic
  !                     sorting.
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
  !                     in the partition of (A,B), LIMIT>=NPTS2
  !                     If LIMIT<NPTS2, the routine will end with
  !                     IER = 6.
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
  !                     IER>0 Abnormal termination of the routine.
  !                             The estimates for integral and error are
  !                             less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more
  !                             subdivisions by increasing the value of
  !                             LIMIT (and taking the according dimension
  !                             adjustments into account). However, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties. If
  !                             the position of a local difficulty can be
  !                             determined (i.e. SINGULARITY,
  !                             DISCONTINUITY within the interval), it
  !                             should be supplied to the routine as an
  !                             element of the vector points. If necessary
  !                             an appropriate special-purpose integrator
  !                             must be used, which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                             The error may be under-estimated.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             At some points of the integration
  !                             interval.
  !                         = 4 The algorithm does not converge.
  !                             Roundoff error is detected in the
  !                             extrapolation table. It is presumed that
  !                             the requested tolerance cannot be
  !                             achieved, and that the returned result is
  !                             the best which can be obtained.
  !                         = 5 The integral is probably divergent, or
  !                             slowly convergent. It must be noted that
  !                             divergence can occur with any other value
  !                             of IER>0.
  !                         = 6 The input is invalid because
  !                             NPTS2<2 or
  !                             Break points are specified outside
  !                             the integration range or
  !                             (EPSABS<=0 and
  !                              EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or LIMIT<NPTS2.
  !                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
  !                             and ELIST(1) are set to zero. ALIST(1) and
  !                             BLIST(1) are set to A and B respectively.
  !
  !            ALIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the left end points
  !                     of the subintervals in the partition of the given
  !                     integration range (A,B)
  !
  !            BLIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the right end points
  !                     of the subintervals in the partition of the given
  !                     integration range (A,B)
  !
  !            RLIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the integral
  !                     approximations on the subintervals
  !
  !            ELIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the moduli of the
  !                     absolute error estimates on the subintervals
  !
  !            PTS    - Double precision
  !                     Vector of dimension at least NPTS2, containing the
  !                     integration limits and the break points of the
  !                     interval in ascending sequence.
  !
  !            LEVEL  - Integer
  !                     Vector of dimension at least LIMIT, containing the
  !                     subdivision levels of the subinterval, i.e. if
  !                     (AA,BB) is a subinterval of (P1,P2) where P1 as
  !                     well as P2 is a user-provided break point or
  !                     integration limit, then (AA,BB) has level L if
  !                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
  !
  !            NDIN   - Integer
  !                     Vector of dimension at least NPTS2, after first
  !                     integration over the intervals (PTS(I)),PTS(I+1),
  !                     I = 0,1, ..., NPTS2-2, the error estimates over
  !                     some of the intervals may have been increased
  !                     artificially, in order to put their subdivision
  !                     forward. If this happens for the subinterval
  !                     numbered K, NDIN(K) is put to 1, otherwise
  !                     NDIN(K) = 0.
  !
  !            IORD   - Integer
  !                     Vector of dimension at least LIMIT, the first K
  !                     elements of which are pointers to the
  !                     error estimates over the subintervals,
  !                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
  !                     form a decreasing sequence, with K = LAST
  !                     If LAST<=(LIMIT/2+2), and K = LIMIT+1-LAST
  !                     otherwise
  !
  !            LAST   - Integer
  !                     Number of subintervals actually produced in the
  !                     subdivisions process
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DQELG, DQK21, DQPSRT

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  INTERFACE
    REAL(DP) FUNCTION F(X)
      IMPORT DP
      REAL(DP) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Ier, Last, Limit, Neval, Npts2
  INTEGER :: Iord(Limit), Level(Limit), Ndin(Npts2)
  REAL(DP) :: A, Abserr, B, Epsabs, Epsrel, Result
  REAL(DP) :: Alist(Limit), Blist(Limit), Elist(Limit), Points(Npts2), Pts(Npts2), &
    Rlist(Limit)
  INTEGER :: i, id, ierro, ind1, ind2, ip1, iroff1, iroff2, iroff3, j, jlow, &
    jupbnd, k, ksgn, ktmin, levcur, levmax, maxerr, nintt, nintp1, npts, &
    nres, nrmax, numrl2
  REAL(DP) :: abseps, area, area1, area12, area2, a1, a2, b1, b2, correc, defabs, &
    defab1, defab2, dres, epmach, erlarg, erlast, errbnd, errmax, error1, erro12, &
    error2, errsum, ertest, oflow, resa, resabs, reseps, res3la(3), rlist2(52), &
    signn, temp, uflow
  LOGICAL :: extrap, noext
  !
  !            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
  !            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
  !            (LIMEXP+2) AT LEAST).
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
  !           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
  !                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
  !                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
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
  !           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
  !                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
  !                       BEEN OBTAINED, IT IS PUT IN RLIST2(NUMRL2) AFTER
  !                       NUMRL2 HAS BEEN INCREASED BY ONE.
  !           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
  !                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
  !           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
  !                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
  !                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
  !                       TRY TO DECREASE THE VALUE OF ERLARG.
  !           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
  !                       NO LONGER ALLOWED (TRUE-VALUE)
  !
  !            MACHINE DEPENDENT CONSTANTS
  !            ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
  !
  !* FIRST EXECUTABLE STATEMENT  DQAGPE
  epmach = D1MACH(4)
  !
  !            TEST ON VALIDITY OF PARAMETERS
  !            -----------------------------
  !
  Ier = 0
  Neval = 0
  Last = 0
  Result = 0._DP
  Abserr = 0._DP
  Alist(1) = A
  Blist(1) = B
  Rlist(1) = 0._DP
  Elist(1) = 0._DP
  Iord(1) = 0
  Level(1) = 0
  npts = Npts2 - 2
  IF( Npts2<2 .OR. Limit<=npts .OR. &
    (Epsabs<=0._DP .AND. Epsrel<MAX(0.5E+02_DP*epmach,0.5E-28_DP)) ) Ier = 6
  IF( Ier==6 ) RETURN
  !
  !            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
  !            ASCENDING SEQUENCE.
  !
  signn = 1._DP
  IF( A>B ) signn = -1._DP
  Pts(1) = MIN(A,B)
  IF( npts/=0 ) THEN
    DO i = 1, npts
      Pts(i+1) = Points(i)
    END DO
  END IF
  Pts(npts+2) = MAX(A,B)
  nintt = npts + 1
  a1 = Pts(1)
  IF( npts/=0 ) THEN
    nintp1 = nintt + 1
    DO i = 1, nintt
      ip1 = i + 1
      DO j = ip1, nintp1
        IF( Pts(i)>Pts(j) ) THEN
          temp = Pts(i)
          Pts(i) = Pts(j)
          Pts(j) = temp
        END IF
      END DO
    END DO
    IF( Pts(1)/=MIN(A,B) .OR. Pts(nintp1)/=MAX(A,B) ) Ier = 6
    IF( Ier==6 ) RETURN
  END IF
  !
  !            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.
  !            ------------------------------------------------
  !
  resabs = 0._DP
  DO i = 1, nintt
    b1 = Pts(i+1)
    CALL DQK21(F,a1,b1,area1,error1,defabs,resa)
    Abserr = Abserr + error1
    Result = Result + area1
    Ndin(i) = 0
    IF( error1==resa .AND. error1/=0._DP ) Ndin(i) = 1
    resabs = resabs + defabs
    Level(i) = 0
    Elist(i) = error1
    Alist(i) = a1
    Blist(i) = b1
    Rlist(i) = area1
    Iord(i) = i
    a1 = b1
  END DO
  errsum = 0._DP
  DO i = 1, nintt
    IF( Ndin(i)==1 ) Elist(i) = Abserr
    errsum = errsum + Elist(i)
  END DO
  !
  !           TEST ON ACCURACY.
  !
  Last = nintt
  Neval = 21*nintt
  dres = ABS(Result)
  errbnd = MAX(Epsabs,Epsrel*dres)
  IF( Abserr<=100._DP*epmach*resabs .AND. Abserr>errbnd ) Ier = 2
  IF( nintt/=1 ) THEN
    DO i = 1, npts
      jlow = i + 1
      ind1 = Iord(i)
      DO j = jlow, nintt
        ind2 = Iord(j)
        IF( Elist(ind1)<=Elist(ind2) ) THEN
          ind1 = ind2
          k = j
        END IF
      END DO
      IF( ind1/=Iord(i) ) THEN
        Iord(k) = Iord(i)
        Iord(i) = ind1
      END IF
    END DO
    IF( Limit<Npts2 ) Ier = 1
  END IF
  IF( Ier/=0 .OR. Abserr<=errbnd ) RETURN
  !
  !           INITIALIZATION
  !           --------------
  !
  rlist2(1) = Result
  maxerr = Iord(1)
  errmax = Elist(maxerr)
  area = Result
  nrmax = 1
  nres = 0
  numrl2 = 1
  ktmin = 0
  extrap = .FALSE.
  noext = .FALSE.
  erlarg = errsum
  ertest = errbnd
  levmax = 1
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ierro = 0
  uflow = D1MACH(1)
  oflow = D1MACH(2)
  Abserr = oflow
  ksgn = -1
  IF( dres>=(1._DP-0.5E+02_DP*epmach)*resabs ) ksgn = 1
  !
  !           MAIN DO-LOOP
  !           ------------
  !
  DO Last = Npts2, Limit
    !
    !           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
    !           ESTIMATE.
    !
    levcur = Level(maxerr) + 1
    a1 = Alist(maxerr)
    b1 = 0.5_DP*(Alist(maxerr)+Blist(maxerr))
    a2 = b1
    b2 = Blist(maxerr)
    erlast = errmax
    CALL DQK21(F,a1,b1,area1,error1,resa,defab1)
    CALL DQK21(F,a2,b2,area2,error2,resa,defab2)
    !
    !           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
    !           AND ERROR AND TEST FOR ACCURACY.
    !
    Neval = Neval + 42
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
    Level(maxerr) = levcur
    Level(Last) = levcur
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
    !           AT A POINT OF THE INTEGRATION RANGE
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
    !- **JUMP OUT OF DO-LOOP
    IF( errsum<=errbnd ) GOTO 200
    !- **JUMP OUT OF DO-LOOP
    IF( Ier/=0 ) EXIT
    IF( .NOT. (noext) ) THEN
      erlarg = erlarg - erlast
      IF( levcur+1<=levmax ) erlarg = erlarg + erro12
      IF( .NOT. (extrap) ) THEN
        !
        !           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
        !           SMALLEST INTERVAL.
        !
        IF( Level(maxerr)+1<=levmax ) CYCLE
        extrap = .TRUE.
        nrmax = 2
      END IF
      IF( ierro/=3 .AND. erlarg>ertest ) THEN
        !
        !           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
        !           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER
        !           THE LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
        !
        id = nrmax
        jupbnd = Last
        IF( Last>(2+Limit/2) ) jupbnd = Limit + 3 - Last
        DO k = id, jupbnd
          maxerr = Iord(nrmax)
          errmax = Elist(maxerr)
          !- **JUMP OUT OF DO-LOOP
          IF( Level(maxerr)+1<=levmax ) GOTO 100
          nrmax = nrmax + 1
        END DO
      END IF
      !
      !           PERFORM EXTRAPOLATION.
      !
      numrl2 = numrl2 + 1
      rlist2(numrl2) = area
      IF( numrl2>2 ) THEN
        CALL DQELG(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin + 1
        IF( ktmin>5 .AND. Abserr<0.1D-02*errsum ) Ier = 5
        IF( abseps<Abserr ) THEN
          ktmin = 0
          Abserr = abseps
          Result = reseps
          correc = erlarg
          ertest = MAX(Epsabs,Epsrel*ABS(reseps))
          !- **JUMP OUT OF DO-LOOP
          IF( Abserr<ertest ) EXIT
        END IF
        !
        !           PREPARE BISECTION OF THE SMALLEST INTERVAL.
        !
        IF( numrl2==1 ) noext = .TRUE.
        IF( Ier>=5 ) EXIT
      END IF
      maxerr = Iord(1)
      errmax = Elist(maxerr)
      nrmax = 1
      extrap = .FALSE.
      levmax = levmax + 1
      erlarg = errsum
    END IF
    100 CONTINUE
  END DO
  !
  !           SET THE FINAL RESULT.
  !           ---------------------
  !
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
    !           TEST ON DIVERGENCE.
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
  300 CONTINUE
  IF( Ier>2 ) Ier = Ier - 1
  Result = Result*signn
  RETURN
END SUBROUTINE DQAGPE
