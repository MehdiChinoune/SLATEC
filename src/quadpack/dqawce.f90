!** DQAWCE
SUBROUTINE DQAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,&
    Alist,Blist,Rlist,Elist,Iord,Last)
  IMPLICIT NONE
  !>
  !***
  !  The routine calculates an approximation result to a
  !            CAUCHY PRINCIPAL VALUE I = Integral of F*W over (A,B)
  !            (W(X) = 1/(X-C), (C.NE.A, C.NE.B), hopefully satisfying
  !            following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I))
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A1, J4
  !***
  ! **Type:**      DOUBLE PRECISION (QAWCE-S, DQAWCE-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE,
  !             CLENSHAW-CURTIS METHOD, QUADPACK, QUADRATURE,
  !             SPECIAL-PURPOSE
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
  !        Computation of a CAUCHY PRINCIPAL VALUE
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
  !            C      - Double precision
  !                     Parameter in the WEIGHT function, C.NE.A, C.NE.B
  !                     If C = A OR C = B, the routine will end with
  !                     IER = 6.
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
  !                     in the partition of (A,B), LIMIT.GE.1
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
  !                             the estimates for integral and error are
  !                             less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more sub-
  !                             divisions by increasing the value of
  !                             LIMIT. However, if this yields no
  !                             improvement it is advised to analyze the
  !                             the integrand, in order to determine the
  !                             the integration difficulties. If the
  !                             position of a local difficulty can be
  !                             determined (e.g. SINGULARITY,
  !                             DISCONTINUITY within the interval) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling
  !                             appropriate integrators on the subranges.
  !                         = 2 The occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 Extremely bad integrand behaviour
  !                             occurs at some interior points of
  !                             the integration interval.
  !                         = 6 The input is invalid, because
  !                             C = A or C = B or
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or LIMIT.LT.1.
  !                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
  !                             IORD(1) and LAST are set to zero. ALIST(1)
  !                             and BLIST(1) are set to A and B
  !                             respectively.
  !
  !            ALIST   - Double precision
  !                      Vector of dimension at least LIMIT, the first
  !                       LAST  elements of which are the left
  !                      end points of the subintervals in the partition
  !                      of the given integration range (A,B)
  !
  !            BLIST   - Double precision
  !                      Vector of dimension at least LIMIT, the first
  !                       LAST  elements of which are the right
  !                      end points of the subintervals in the partition
  !                      of the given integration range (A,B)
  !
  !            RLIST   - Double precision
  !                      Vector of dimension at least LIMIT, the first
  !                       LAST  elements of which are the integral
  !                      approximations on the subintervals
  !
  !            ELIST   - Double precision
  !                      Vector of dimension LIMIT, the first  LAST
  !                      elements of which are the moduli of the absolute
  !                      error estimates on the subintervals
  !
  !            IORD    - Integer
  !                      Vector of dimension at least LIMIT, the first K
  !                      elements of which are pointers to the error
  !                      estimates over the subintervals, so that
  !                      ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
  !                      If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
  !                      otherwise, form a decreasing sequence
  !
  !            LAST    - Integer
  !                      Number of subintervals actually produced in
  !                      the subdivision process
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DQC25C, DQPSRT

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  !
  REAL(8) :: A, aa, Abserr, Alist(*), area, area1, area12, area2, &
    a1, a2, B, bb, Blist(*), b1, b2, C, Elist(*), &
    epmach, Epsabs, Epsrel, errbnd, errmax, error1, &
    erro12, error2, errsum, Result, Rlist(*), uflow
  INTEGER Ier, Iord(*), iroff1, iroff2, k, krule, Last, Limit, maxerr, &
    nev, Neval, nrmax
  !
  REAL(8), EXTERNAL :: F
  REAL(8), EXTERNAL :: D1MACH
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
  !* FIRST EXECUTABLE STATEMENT  DQAWCE
  epmach = D1MACH(4)
  uflow = D1MACH(1)
  !
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  Ier = 6
  Neval = 0
  Last = 0
  Alist(1) = A
  Blist(1) = B
  Rlist(1) = 0.0D+00
  Elist(1) = 0.0D+00
  Iord(1) = 0
  Result = 0.0D+00
  Abserr = 0.0D+00
  IF ( .NOT.(C==A.OR.C==B.OR.(Epsabs<=0.0D+00.AND.Epsrel<MAX(0.5D+02*epmach,&
      0.5D-28))) ) THEN
    !
    !           FIRST APPROXIMATION TO THE INTEGRAL
    !           -----------------------------------
    !
    aa = A
    bb = B
    IF ( A>B ) THEN
      aa = B
      bb = A
    ENDIF
    Ier = 0
    krule = 1
    CALL DQC25C(F,aa,bb,C,Result,Abserr,krule,Neval)
    Last = 1
    Rlist(1) = Result
    Elist(1) = Abserr
    Iord(1) = 1
    Alist(1) = A
    Blist(1) = B
    !
    !           TEST ON ACCURACY
    !
    errbnd = MAX(Epsabs,Epsrel*ABS(Result))
    IF ( Limit==1 ) Ier = 1
    IF ( Abserr>=MIN(0.1D-01*ABS(Result),errbnd).AND.Ier/=1 ) THEN
      !
      !           INITIALIZATION
      !           --------------
      !
      Alist(1) = aa
      Blist(1) = bb
      Rlist(1) = Result
      errmax = Abserr
      maxerr = 1
      area = Result
      errsum = Abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
      !
      !           MAIN DO-LOOP
      !           ------------
      !
      DO Last = 2, Limit
        !
        !           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST
        !           ERROR ESTIMATE.
        !
        a1 = Alist(maxerr)
        b1 = 0.5D+00*(Alist(maxerr)+Blist(maxerr))
        b2 = Blist(maxerr)
        IF ( C<=b1.AND.C>a1 ) b1 = 0.5D+00*(C+b2)
        IF ( C>b1.AND.C<b2 ) b1 = 0.5D+00*(a1+C)
        a2 = b1
        krule = 2
        CALL DQC25C(F,a1,b1,C,area1,error1,krule,nev)
        Neval = Neval + nev
        CALL DQC25C(F,a2,b2,C,area2,error2,krule,nev)
        Neval = Neval + nev
        !
        !           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
        !           AND ERROR AND TEST FOR ACCURACY.
        !
        area12 = area1 + area2
        erro12 = error1 + error2
        errsum = errsum + erro12 - errmax
        area = area + area12 - Rlist(maxerr)
        IF ( ABS(Rlist(maxerr)-area12)<0.1D-04*ABS(area12).AND.&
          erro12>=0.99D+00*errmax.AND.krule==0 ) iroff1 = iroff1 + 1
        IF ( Last>10.AND.erro12>errmax.AND.krule==0 ) iroff2 = iroff2 + 1
        Rlist(maxerr) = area1
        Rlist(Last) = area2
        errbnd = MAX(Epsabs,Epsrel*ABS(area))
        IF ( errsum>errbnd ) THEN
          !
          !           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
          !
          IF ( iroff1>=6.AND.iroff2>20 ) Ier = 2
          !
          !           SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL
          !           BISECTIONS EXCEEDS LIMIT.
          !
          IF ( Last==Limit ) Ier = 1
          !
          !           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
          !           AT A POINT OF THE INTEGRATION RANGE.
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
        !           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
        !
        CALL DQPSRT(Limit,Last,maxerr,errmax,Elist,Iord,nrmax)
        !- **JUMP OUT OF DO-LOOP
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
    IF ( aa==B ) Result = -Result
  ENDIF
END SUBROUTINE DQAWCE
