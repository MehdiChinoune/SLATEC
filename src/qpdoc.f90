!DECK QPDOC
SUBROUTINE QPDOC
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QPDOC
  !***PURPOSE  Documentation for QUADPACK, a package of subprograms for
  !            automatic evaluation of one-dimensional definite integrals.
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2, Z
  !***TYPE      ALL (QPDOC-A)
  !***KEYWORDS  DOCUMENTATION, GUIDELINES FOR SELECTION, QUADPACK,
  !             QUADRATURE, SURVEY OF INTEGRATORS
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           Kahaner, D. K., (NBS)
  !***DESCRIPTION
  !
  ! 1. Introduction
  !    ------------
  !    QUADPACK is a FORTRAN subroutine package for the numerical
  !    computation of definite one-dimensional integrals. It originated
  !    from a joint project of R. Piessens and E. de Doncker (Appl.
  !    Math. and Progr. Div.- K.U.Leuven, Belgium), C. Ueberhuber (Inst.
  !    Fuer Math.- Techn. U. Wien, Austria), and D. Kahaner (National
  !    Bureau of Standards- Washington D.C., U.S.A.).
  !
  !    Documentation routine QPDOC describes the package in the form it
  !    was released from A.M.P.D.- Leuven, for adherence to the SLATEC
  !    library in May 1981. Apart from a survey of the integrators, some
  !    guidelines will be given in order to help the QUADPACK user with
  !    selecting an appropriate routine or a combination of several
  !    routines for handling his problem.
  !
  !    In the Long Description of QPDOC it is demonstrated how to call
  !    the integrators, by means of small example calling programs.
  !
  !    For precise guidelines involving the use of each routine in
  !    particular, we refer to the extensive introductory comments
  !    within each routine.
  !
  ! 2. Survey
  !    ------
  !    The following list gives an overview of the QUADPACK integrators.
  !    The routine names for the DOUBLE PRECISION versions are preceded
  !    by the letter D.
  !
  !    - QNG  : Is a simple non-adaptive automatic integrator, based on
  !             a sequence of rules with increasing degree of algebraic
  !             precision (Patterson, 1968).
  !
  !    - QAG  : Is a simple globally adaptive integrator using the
  !             strategy of Aind (Piessens, 1973). It is possible to
  !             choose between 6 pairs of Gauss-Kronrod quadrature
  !             formulae for the rule evaluation component. The pairs
  !             of high degree of precision are suitable for handling
  !             integration difficulties due to a strongly oscillating
  !             integrand.
  !
  !    - QAGS : Is an integrator based on globally adaptive interval
  !             subdivision in connection with extrapolation (de Doncker,
  !             1978) by the Epsilon algorithm (Wynn, 1956).
  !
  !    - QAGP : Serves the same purposes as QAGS, but also allows
  !             for eventual user-supplied information, i.e. the
  !             abscissae of internal singularities, discontinuities
  !             and other difficulties of the integrand function.
  !             The algorithm is a modification of that in QAGS.
  !
  !    - QAGI : Handles integration over infinite intervals. The
  !             infinite range is mapped onto a finite interval and
  !             then the same strategy as in QAGS is applied.
  !
  !    - QAWO : Is a routine for the integration of COS(OMEGA*X)*F(X)
  !             or SIN(OMEGA*X)*F(X) over a finite interval (A,B).
  !             OMEGA is is specified by the user
  !             The rule evaluation component is based on the
  !             modified Clenshaw-Curtis technique.
  !             An adaptive subdivision scheme is used connected with
  !             an extrapolation procedure, which is a modification
  !             of that in QAGS and provides the possibility to deal
  !             even with singularities in F.
  !
  !    - QAWF : Calculates the Fourier cosine or Fourier sine
  !             transform of F(X), for user-supplied interval (A,
  !             INFINITY), OMEGA, and F. The procedure of QAWO is
  !             used on successive finite intervals, and convergence
  !             acceleration by means of the Epsilon algorithm (Wynn,
  !             1956) is applied to the series of the integral
  !             contributions.
  !
  !    - QAWS : Integrates W(X)*F(X) over (A,B) with A.LT.B finite,
  !             and   W(X) = ((X-A)**ALFA)*((B-X)**BETA)*V(X)
  !             where V(X) = 1 or LOG(X-A) or LOG(B-X)
  !                            or LOG(X-A)*LOG(B-X)
  !             and   ALFA.GT.(-1), BETA.GT.(-1).
  !             The user specifies A, B, ALFA, BETA and the type of
  !             the function V.
  !             A globally adaptive subdivision strategy is applied,
  !             with modified Clenshaw-Curtis integration on the
  !             subintervals which contain A or B.
  !
  !    - QAWC : Computes the Cauchy Principal Value of F(X)/(X-C)
  !             over a finite interval (A,B) and for
  !             user-determined C.
  !             The strategy is globally adaptive, and modified
  !             Clenshaw-Curtis integration is used on the subranges
  !             which contain the point X = C.
  !
  !  Each of the routines above also has a "more detailed" version
  !    with a name ending in E, as QAGE.  These provide more
  !    information and control than the easier versions.
  !
  !
  !   The preceding routines are all automatic.  That is, the user
  !      inputs his problem and an error tolerance.  The routine
  !      attempts to perform the integration to within the requested
  !      absolute or relative error.
  !   There are, in addition, a number of non-automatic integrators.
  !      These are most useful when the problem is such that the
  !      user knows that a fixed rule will provide the accuracy
  !      required.  Typically they return an error estimate but make
  !      no attempt to satisfy any particular input error request.
  !
  !      QK15
  !      QK21
  !      QK31
  !      QK41
  !      QK51
  !      QK61
  !           Estimate the integral on [a,b] using 15, 21,..., 61
  !           point rule and return an error estimate.
  !      QK15I 15 point rule for (semi)infinite interval.
  !      QK15W 15 point rule for special singular weight functions.
  !      QC25C 25 point rule for Cauchy Principal Values
  !      QC25F 25 point rule for sin/cos integrand.
  !      QMOMO Integrates k-th degree Chebyshev polynomial times
  !            function with various explicit singularities.
  !
  ! 3. Guidelines for the use of QUADPACK
  !    ----------------------------------
  !    Here it is not our purpose to investigate the question when
  !    automatic quadrature should be used. We shall rather attempt
  !    to help the user who already made the decision to use QUADPACK,
  !    with selecting an appropriate routine or a combination of
  !    several routines for handling his problem.
  !
  !    For both quadrature over finite and over infinite intervals,
  !    one of the first questions to be answered by the user is
  !    related to the amount of computer time he wants to spend,
  !    versus his -own- time which would be needed, for example, for
  !    manual subdivision of the interval or other analytic
  !    manipulations.
  !
  !    (1) The user may not care about computer time, or not be
  !        willing to do any analysis of the problem. especially when
  !        only one or a few integrals must be calculated, this attitude
  !        can be perfectly reasonable. In this case it is clear that
  !        either the most sophisticated of the routines for finite
  !        intervals, QAGS, must be used, or its analogue for infinite
  !        intervals, GAGI. These routines are able to cope with
  !        rather difficult, even with improper integrals.
  !        This way of proceeding may be expensive. But the integrator
  !        is supposed to give you an answer in return, with additional
  !        information in the case of a failure, through its error
  !        estimate and flag. Yet it must be stressed that the programs
  !        cannot be totally reliable.
  !        ------
  !
  !    (2) The user may want to examine the integrand function.
  !        If bad local difficulties occur, such as a discontinuity, a
  !        singularity, derivative singularity or high peak at one or
  !        more points within the interval, the first advice is to
  !        split up the interval at these points. The integrand must
  !        then be examined over each of the subintervals separately,
  !        so that a suitable integrator can be selected for each of
  !        them. If this yields problems involving relative accuracies
  !        to be imposed on -finite- subintervals, one can make use of
  !        QAGP, which must be provided with the positions of the local
  !        difficulties. However, if strong singularities are present
  !        and a high accuracy is requested, application of QAGS on the
  !        subintervals may yield a better result.
  !
  !        For quadrature over finite intervals we thus dispose of QAGS
  !        and
  !        - QNG for well-behaved integrands,
  !        - QAG for functions with an oscillating behaviour of a non
  !          specific type,
  !        - QAWO for functions, eventually singular, containing a
  !          factor COS(OMEGA*X) or SIN(OMEGA*X) where OMEGA is known,
  !        - QAWS for integrands with Algebraico-Logarithmic end point
  !          singularities of known type,
  !        - QAWC for Cauchy Principal Values.
  !
  !        Remark
  !        ------
  !        On return, the work arrays in the argument lists of the
  !        adaptive integrators contain information about the interval
  !        subdivision process and hence about the integrand behaviour:
  !        the end points of the subintervals, the local integral
  !        contributions and error estimates, and eventually other
  !        characteristics. For this reason, and because of its simple
  !        globally adaptive nature, the routine QAG in particular is
  !        well-suited for integrand examination. Difficult spots can
  !        be located by investigating the error estimates on the
  !        subintervals.
  !
  !        For infinite intervals we provide only one general-purpose
  !        routine, QAGI. It is based on the QAGS algorithm applied
  !        after a transformation of the original interval into (0,1).
  !        Yet it may eventuate that another type of transformation is
  !        more appropriate, or one might prefer to break up the
  !        original interval and use QAGI only on the infinite part
  !        and so on. These kinds of actions suggest a combined use of
  !        different QUADPACK integrators. Note that, when the only
  !        difficulty is an integrand singularity at the finite
  !        integration limit, it will in general not be necessary to
  !        break up the interval, as QAGI deals with several types of
  !        singularity at the boundary point of the integration range.
  !        It also handles slowly convergent improper integrals, on
  !        the condition that the integrand does not oscillate over
  !        the entire infinite interval. If it does we would advise
  !        to sum succeeding positive and negative contributions to
  !        the integral -e.g. integrate between the zeros- with one
  !        or more of the finite-range integrators, and apply
  !        convergence acceleration eventually by means of QUADPACK
  !        subroutine QELG which implements the Epsilon algorithm.
  !        Such quadrature problems include the Fourier transform as
  !        a special case. Yet for the latter we have an automatic
  !        integrator available, QAWF.
  !
  ! *Long Description:
  !
  ! 4. Example Programs
  !    ----------------
  ! 4.1. Calling Program for QNG
  !      -----------------------
  !
  !            REAL A,ABSERR,B,F,EPSABS,EPSREL,RESULT
  !            INTEGER IER,NEVAL
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            CALL QNG(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = EXP(X)/(X*X+0.1E+01)
  !            RETURN
  !            END
  !
  ! 4.2. Calling Program for QAG
  !      -----------------------
  !
  !            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
  !            INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
  !            DIMENSION IWORK(100),WORK(400)
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            KEY = 6
  !            LIMIT = 100
  !            LENW = LIMIT*4
  !            CALL QAG(F,A,B,EPSABS,EPSREL,KEY,RESULT,ABSERR,NEVAL,
  !           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 2.0E0/(2.0E0+SIN(31.41592653589793E0*X))
  !            RETURN
  !            END
  !
  ! 4.3. Calling Program for QAGS
  !      ------------------------
  !
  !            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
  !            INTEGER IER,IWORK,LAST,LENW,LIMIT,NEVAL
  !            DIMENSION IWORK(100),WORK(400)
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            LIMIT = 100
  !            LENW = LIMIT*4
  !            CALL QAGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,IER,
  !           *  LIMIT,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 0.0E0
  !            IF(X.GT.0.0E0) F = 1.0E0/SQRT(X)
  !            RETURN
  !            END
  !
  ! 4.4. Calling Program for QAGP
  !      ------------------------
  !
  !            REAL A,ABSERR,B,EPSABS,EPSREL,F,POINTS,RESULT,WORK
  !            INTEGER IER,IWORK,LAST,LENIW,LENW,LIMIT,NEVAL,NPTS2
  !            DIMENSION IWORK(204),POINTS(4),WORK(404)
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            NPTS2 = 4
  !            POINTS(1) = 1.0E0/7.0E0
  !            POINTS(2) = 2.0E0/3.0E0
  !            LIMIT = 100
  !            LENIW = LIMIT*2+NPTS2
  !            LENW = LIMIT*4+NPTS2
  !            CALL QAGP(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,RESULT,ABSERR,
  !           *  NEVAL,IER,LENIW,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 0.0E+00
  !            IF(X.NE.1.0E0/7.0E0.AND.X.NE.2.0E0/3.0E0) F =
  !           *  ABS(X-1.0E0/7.0E0)**(-0.25E0)*
  !           *  ABS(X-2.0E0/3.0E0)**(-0.55E0)
  !            RETURN
  !            END
  !
  ! 4.5. Calling Program for QAGI
  !      ------------------------
  !
  !            REAL ABSERR,BOUN,EPSABS,EPSREL,F,RESULT,WORK
  !            INTEGER IER,INF,IWORK,LAST,LENW,LIMIT,NEVAL
  !            DIMENSION IWORK(100),WORK(400)
  !            EXTERNAL F
  !            BOUN = 0.0E0
  !            INF = 1
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            LIMIT = 100
  !            LENW = LIMIT*4
  !            CALL QAGI(F,BOUN,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
  !           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 0.0E0
  !            IF(X.GT.0.0E0) F = SQRT(X)*LOG(X)/
  !           *                   ((X+1.0E0)*(X+2.0E0))
  !            RETURN
  !            END
  !
  ! 4.6. Calling Program for QAWO
  !      ------------------------
  !
  !            REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,OMEGA,WORK
  !            INTEGER IER,INTEGR,IWORK,LAST,LENIW,LENW,LIMIT,MAXP1,NEVAL
  !            DIMENSION IWORK(200),WORK(925)
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            OMEGA = 10.0E0
  !            INTEGR = 1
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            LIMIT = 100
  !            LENIW = LIMIT*2
  !            MAXP1 = 21
  !            LENW = LIMIT*4+MAXP1*25
  !            CALL QAWO(F,A,B,OMEGA,INTEGR,EPSABS,EPSREL,RESULT,ABSERR,
  !           *  NEVAL,IER,LENIW,MAXP1,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 0.0E0
  !            IF(X.GT.0.0E0) F = EXP(-X)*LOG(X)
  !            RETURN
  !            END
  !
  ! 4.7. Calling Program for QAWF
  !      ------------------------
  !
  !            REAL A,ABSERR,EPSABS,F,RESULT,OMEGA,WORK
  !            INTEGER IER,INTEGR,IWORK,LAST,LENIW,LENW,LIMIT,LIMLST,
  !           *  LST,MAXP1,NEVAL
  !            DIMENSION IWORK(250),WORK(1025)
  !            EXTERNAL F
  !            A = 0.0E0
  !            OMEGA = 8.0E0
  !            INTEGR = 2
  !            EPSABS = 1.0E-3
  !            LIMLST = 50
  !            LIMIT = 100
  !            LENIW = LIMIT*2+LIMLST
  !            MAXP1 = 21
  !            LENW = LENIW*2+MAXP1*25
  !            CALL QAWF(F,A,OMEGA,INTEGR,EPSABS,RESULT,ABSERR,NEVAL,
  !           *  IER,LIMLST,LST,LENIW,MAXP1,LENW,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            IF(X.GT.0.0E0) F = SIN(50.0E0*X)/(X*SQRT(X))
  !            RETURN
  !            END
  !
  ! 4.8. Calling Program for QAWS
  !      ------------------------
  !
  !            REAL A,ABSERR,ALFA,B,BETA,EPSABS,EPSREL,F,RESULT,WORK
  !            INTEGER IER,INTEGR,IWORK,LAST,LENW,LIMIT,NEVAL
  !            DIMENSION IWORK(100),WORK(400)
  !            EXTERNAL F
  !            A = 0.0E0
  !            B = 1.0E0
  !            ALFA = -0.5E0
  !            BETA = -0.5E0
  !            INTEGR = 1
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            LIMIT = 100
  !            LENW = LIMIT*4
  !            CALL QAWS(F,A,B,ALFA,BETA,INTEGR,EPSABS,EPSREL,RESULT,
  !           *  ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = SIN(10.0E0*X)
  !            RETURN
  !            END
  !
  ! 4.9. Calling Program for QAWC
  !      ------------------------
  !
  !            REAL A,ABSERR,B,C,EPSABS,EPSREL,F,RESULT,WORK
  !            INTEGER IER,IWORK,LAST,LENW,LIMIT,NEVAL
  !            DIMENSION IWORK(100),WORK(400)
  !            EXTERNAL F
  !            A = -1.0E0
  !            B = 1.0E0
  !            C = 0.5E0
  !            EPSABS = 0.0E0
  !            EPSREL = 1.0E-3
  !            LIMIT = 100
  !            LENW = LIMIT*4
  !            CALL QAWC(F,A,B,C,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
  !           *  IER,LIMIT,LENW,LAST,IWORK,WORK)
  !      C  INCLUDE WRITE STATEMENTS
  !            STOP
  !            END
  !      C
  !            REAL FUNCTION F(X)
  !            REAL X
  !            F = 1.0E0/(X*X+1.0E-4)
  !            RETURN
  !            END
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900723  PURPOSE section revised.  (WRB)
  !***END PROLOGUE  QPDOC
  !***FIRST EXECUTABLE STATEMENT  QPDOC
END SUBROUTINE QPDOC
