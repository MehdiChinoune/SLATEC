!*==DQAWSE.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQAWSE
SUBROUTINE DQAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit,Result,&
    Abserr,Neval,Ier,Alist,Blist,Rlist,Elist,Iord,Last)
  IMPLICIT NONE
  !*--DQAWSE6
  !***BEGIN PROLOGUE  DQAWSE
  !***PURPOSE  The routine calculates an approximation result to a given
  !            definite integral I = Integral of F*W over (A,B),
  !            (where W shows a singular behaviour at the end points,
  !            see parameter INTEGR).
  !            Hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A2A1
  !***TYPE      DOUBLE PRECISION (QAWSE-S, DQAWSE-D)
  !***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES,
  !             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, QUADPACK,
  !             QUADRATURE, SPECIAL-PURPOSE
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !        Integration of functions having algebraico-logarithmic
  !        end point singularities
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
  !                     Upper limit of integration, B.GT.A
  !                     If B.LE.A, the routine will end with IER = 6.
  !
  !            ALFA   - Double precision
  !                     Parameter in the WEIGHT function, ALFA.GT.(-1)
  !                     If ALFA.LE.(-1), the routine will end with
  !                     IER = 6.
  !
  !            BETA   - Double precision
  !                     Parameter in the WEIGHT function, BETA.GT.(-1)
  !                     If BETA.LE.(-1), the routine will end with
  !                     IER = 6.
  !
  !            INTEGR - Integer
  !                     Indicates which WEIGHT function is to be used
  !                     = 1  (X-A)**ALFA*(B-X)**BETA
  !                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
  !                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
  !                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
  !                     If INTEGR.LT.1 or INTEGR.GT.4, the routine
  !                     will end with IER = 6.
  !
  !            EPSABS - Double precision
  !                     Absolute accuracy requested
  !            EPSREL - Double precision
  !                     Relative accuracy requested
  !                     If  EPSABS.LE.0
  !                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     the routine will end with IER = 6.
  !
  !            LIMIT  - Integer
  !                     Gives an upper bound on the number of subintervals
  !                     in the partition of (A,B), LIMIT.GE.2
  !                     If LIMIT.LT.2, the routine will end with IER = 6.
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
  !                     IER.GT.0 Abnormal termination of the routine
  !                             the estimates for the integral and error
  !                             are less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                         = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more
  !                             subdivisions by increasing the value of
  !                             LIMIT. However, if this yields no
  !                             improvement, it is advised to analyze the
  !                             integrand in order to determine the
  !                             integration difficulties which prevent the
  !                             requested tolerance from being achieved.
  !                             In case of a jump DISCONTINUITY or a local
  !                             SINGULARITY of algebraico-logarithmic type
  !                             at one or more interior points of the
  !                             integration range, one should proceed by
  !                             splitting up the interval at these
  !                             points and calling the integrator on the
  !                             subranges.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 The input is invalid, because
  !                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1), or
  !                             INTEGR.LT.1 or INTEGR.GT.4, or
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                             or LIMIT.LT.2.
  !                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
  !                             IORD(1) and LAST are set to zero. ALIST(1)
  !                             and BLIST(1) are set to A and B
  !                             respectively.
  !
  !            ALIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the left
  !                     end points of the subintervals in the partition
  !                     of the given integration range (A,B)
  !
  !            BLIST  - Double precision
  !                     Vector of dimension at least LIMIT, the first
  !                      LAST  elements of which are the right
  !                     end points of the subintervals in the partition
  !                     of the given integration range (A,B)
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
  !            IORD   - Integer
  !                     Vector of dimension at least LIMIT, the first K
  !                     of which are pointers to the error
  !                     estimates over the subintervals, so that
  !                     ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
  !                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
  !                     otherwise form a decreasing sequence
  !
  !            LAST   - Integer
  !                     Number of subintervals actually produced in
  !                     the subdivision process
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DQC25S, DQMOMO, DQPSRT
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DQAWSE
  !
  REAL(8) :: A, Abserr, Alfa, Alist, area, area1, area12, &
    area2, a1, a2, B, Beta, Blist, b1, b2, centre, &
    D1MACH, Elist, epmach, Epsabs, Epsrel, errbnd, &
    errmax, error1, erro12, error2, errsum, F, resas1, &
    resas2, Result, rg, rh, ri, rj, Rlist, uflow
  INTEGER Ier, Integr, Iord, iroff1, iroff2, k, Last, Limit, &
    maxerr, nev, Neval, nrmax
  !
  EXTERNAL F
  !
  DIMENSION Alist(*), Blist(*), Rlist(*), Elist(*), Iord(*), ri(25), &
    rj(25), rh(25), rg(25)
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
  !           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
  !           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
  !                       ERROR ESTIMATE
  !           ERRMAX    - ELIST(MAXERR)
  !           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
  !           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
  !           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
  !                       ABS(RESULT))
  !           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
  !           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
  !           LAST      - INDEX FOR SUBDIVISION
  !
  !
  !            MACHINE DEPENDENT CONSTANTS
  !            ---------------------------
  !
  !           EPMACH IS THE LARGEST RELATIVE SPACING.
  !           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
  !
  !***FIRST EXECUTABLE STATEMENT  DQAWSE
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  Ier = 6
  Neval = 0
  Last = 0
  Rlist(1) = 0.0D+00
  Elist(1) = 0.0D+00
  Iord(1) = 0
  Result = 0.0D+00
  Abserr = 0.0D+00
  IF ( .NOT.(B<=A.OR.(Epsabs==0.0D+00.AND.Epsrel<MAX(0.5D+02*epmach,0.5D-28)&
      ).OR.Alfa<=(-0.1D+01).OR.Beta<=(-0.1D+01).OR.Integr<1.OR.Integr>4.OR.&
      Limit<2) ) THEN
    Ier = 0
    !
    !           COMPUTE THE MODIFIED CHEBYSHEV MOMENTS.
    !
    CALL DQMOMO(Alfa,Beta,ri,rj,rg,rh,Integr)
    !
    !           INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B).
    !
    centre = 0.5D+00*(B+A)
    CALL DQC25S(F,A,B,A,centre,Alfa,Beta,ri,rj,rg,rh,area1,error1,resas1,&
      Integr,nev)
    Neval = nev
    CALL DQC25S(F,A,B,centre,B,Alfa,Beta,ri,rj,rg,rh,area2,error2,resas2,&
      Integr,nev)
    Last = 2
    Neval = Neval + nev
    Result = area1 + area2
    Abserr = error1 + error2
    !
    !           TEST ON ACCURACY.
    !
    errbnd = MAX(Epsabs,Epsrel*ABS(Result))
    !
    !           INITIALIZATION
    !           --------------
    !
    IF ( error2>error1 ) THEN
      Alist(1) = centre
      Alist(2) = A
      Blist(1) = B
      Blist(2) = centre
      Rlist(1) = area2
      Rlist(2) = area1
      Elist(1) = error2
      Elist(2) = error1
    ELSE
      Alist(1) = A
      Alist(2) = centre
      Blist(1) = centre
      Blist(2) = B
      Rlist(1) = area1
      Rlist(2) = area2
      Elist(1) = error1
      Elist(2) = error2
    ENDIF
    Iord(1) = 1
    Iord(2) = 2
    IF ( Limit==2 ) Ier = 1
    IF ( Abserr>errbnd.AND.Ier/=1 ) THEN
      errmax = Elist(1)
      maxerr = 1
      nrmax = 1
      area = Result
      errsum = Abserr
      iroff1 = 0
      iroff2 = 0
      !
      !            MAIN DO-LOOP
      !            ------------
      !
      DO Last = 3, Limit
        !
        !           BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE.
        !
        a1 = Alist(maxerr)
        b1 = 0.5D+00*(Alist(maxerr)+Blist(maxerr))
        a2 = b1
        b2 = Blist(maxerr)
        !
        CALL DQC25S(F,A,B,a1,b1,Alfa,Beta,ri,rj,rg,rh,area1,error1,resas1,&
          Integr,nev)
        Neval = Neval + nev
        CALL DQC25S(F,A,B,a2,b2,Alfa,Beta,ri,rj,rg,rh,area2,error2,resas2,&
          Integr,nev)
        Neval = Neval + nev
        !
        !           IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR
        !           AND TEST FOR ACCURACY.
        !
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 - errmax
        area = area + area12 - Rlist(maxerr)
        IF ( A/=a1.AND.B/=b2 ) THEN
          IF ( resas1/=error1.AND.resas2/=error2 ) THEN
            !
            !           TEST FOR ROUNDOFF ERROR.
            !
            IF ( ABS(Rlist(maxerr)-area12)<0.1D-04*ABS(area12).AND.&
              erro12>=0.99D+00*errmax ) iroff1 = iroff1 + 1
            IF ( Last>10.AND.erro12>errmax ) iroff2 = iroff2 + 1
          ENDIF
        ENDIF
        Rlist(maxerr) = area1
        Rlist(Last) = area2
        !
        !           TEST ON ACCURACY.
        !
        errbnd = MAX(Epsabs,Epsrel*ABS(area))
        IF ( errsum>errbnd ) THEN
          !
          !           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL
          !           BISECTIONS EXCEEDS LIMIT.
          !
          IF ( Last==Limit ) Ier = 1
          !
          !
          !           SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR.
          !
          IF ( iroff1>=6.OR.iroff2>=20 ) Ier = 2
          !
          !           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
          !           AT INTERIOR POINTS OF INTEGRATION RANGE.
          !
          IF ( MAX(ABS(a1),ABS(b2))<=(0.1D+01+0.1D+03*epmach)&
            *(ABS(a2)+0.1D+04*uflow) ) Ier = 3
        ENDIF
        !
        !           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
        !
        IF ( error2>error1 ) THEN
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
        ENDIF
        !
        !           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
        !           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
        !           WITH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
        !
        CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
        ! ***JUMP OUT OF DO-LOOP
        IF ( Ier/=0.OR.errsum<=errbnd ) EXIT
      ENDDO
      !
      !           COMPUTE FINAL RESULT.
      !           ---------------------
      !
      Result = 0.0D+00
      DO k = 1, Last
        Result = Result + Rlist(k)
      ENDDO
      Abserr = errsum
    ENDIF
  ENDIF
END SUBROUTINE DQAWSE
