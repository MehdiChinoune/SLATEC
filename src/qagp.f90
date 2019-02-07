!*==QAGP.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QAGP
SUBROUTINE QAGP(F,A,B,Npts2,Points,Epsabs,Epsrel,Result,Abserr,Neval,Ier,&
    Leniw,Lenw,Last,Iwork,Work)
  IMPLICIT NONE
  !*--QAGP6
  !*** Start of declarations inserted by SPAG
  INTEGER l4 , Last
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  QAGP
  !***PURPOSE  The routine calculates an approximation result to a given
  !            definite integral I = Integral of F over (A,B),
  !            hopefully satisfying following claim for accuracy
  !            break points of the integration interval, where local
  !            difficulties of the integrand may occur(e.g. SINGULARITIES,
  !            DISCONTINUITIES), are provided by the user.
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A2A1
  !***TYPE      SINGLE PRECISION (QAGP-S, DQAGP-D)
  !***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
  !             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
  !             SINGULARITIES AT USER SPECIFIED POINTS
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !        Computation of a definite integral
  !        Standard fortran subroutine
  !        Real version
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Real
  !                     Function subprogram defining the integrand
  !                     Function F(X). The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            A      - Real
  !                     Lower limit of integration
  !
  !            B      - Real
  !                     Upper limit of integration
  !
  !            NPTS2  - Integer
  !                     Number equal to two more than the number of
  !                     user-supplied break points within the integration
  !                     range, NPTS.GE.2.
  !                     If NPTS2.LT.2, The routine will end with IER = 6.
  !
  !            POINTS - Real
  !                     Vector of dimension NPTS2, the first (NPTS2-2)
  !                     elements of which are the user provided break
  !                     points. If these points do not constitute an
  !                     ascending sequence there will be an automatic
  !                     sorting.
  !
  !            EPSABS - Real
  !                     Absolute accuracy requested
  !            EPSREL - Real
  !                     Relative accuracy requested
  !                     If  EPSABS.LE.0
  !                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     The routine will end with IER = 6.
  !
  !         ON RETURN
  !            RESULT - Real
  !                     Approximation to the integral
  !
  !            ABSERR - Real
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
  !                     IER.GT.0 Abnormal termination of the routine.
  !                             The estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. one can allow more
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
  !                             at some points of the integration
  !                             interval.
  !                         = 4 The algorithm does not converge.
  !                             roundoff error is detected in the
  !                             extrapolation table.
  !                             It is presumed that the requested
  !                             tolerance cannot be achieved, and that
  !                             the returned RESULT is the best which
  !                             can be obtained.
  !                         = 5 The integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of IER.GT.0.
  !                         = 6 The input is invalid because
  !                             NPTS2.LT.2 or
  !                             break points are specified outside
  !                             the integration range or
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             RESULT, ABSERR, NEVAL, LAST are set to
  !                             zero.  Except when LENIW or LENW or NPTS2
  !                             is invalid, IWORK(1), IWORK(LIMIT+1),
  !                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1)
  !                             are set to zero.
  !                             WORK(1) is set to A and WORK(LIMIT+1)
  !                             to B (where LIMIT = (LENIW-NPTS2)/2).
  !
  !         DIMENSIONING PARAMETERS
  !            LENIW - Integer
  !                    Dimensioning parameter for IWORK
  !                    LENIW determines LIMIT = (LENIW-NPTS2)/2,
  !                    which is the maximum number of subintervals in the
  !                    partition of the given integration interval (A,B),
  !                    LENIW.GE.(3*NPTS2-2).
  !                    If LENIW.LT.(3*NPTS2-2), the routine will end with
  !                    IER = 6.
  !
  !            LENW  - Integer
  !                    Dimensioning parameter for WORK
  !                    LENW must be at least LENIW*2-NPTS2.
  !                    If LENW.LT.LENIW*2-NPTS2, the routine will end
  !                    with IER = 6.
  !
  !            LAST  - Integer
  !                    On return, LAST equals the number of subintervals
  !                    produced in the subdivision process, which
  !                    determines the number of significant elements
  !                    actually in the WORK ARRAYS.
  !
  !         WORK ARRAYS
  !            IWORK - Integer
  !                    Vector of dimension at least LENIW. on return,
  !                    the first K elements of which contain
  !                    pointers to the error estimates over the
  !                    subintervals, such that WORK(LIMIT*3+IWORK(1)),...,
  !                    WORK(LIMIT*3+IWORK(K)) form a decreasing
  !                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
  !                    K = LIMIT+1-LAST otherwise
  !                    IWORK(LIMIT+1), ...,IWORK(LIMIT+LAST) Contain the
  !                     subdivision levels of the subintervals, i.e.
  !                     if (AA,BB) is a subinterval of (P1,P2)
  !                     where P1 as well as P2 is a user-provided
  !                     break point or integration LIMIT, then (AA,BB) has
  !                     level L if ABS(BB-AA) = ABS(P2-P1)*2**(-L),
  !                    IWORK(LIMIT*2+1), ..., IWORK(LIMIT*2+NPTS2) have
  !                     no significance for the user,
  !                    note that LIMIT = (LENIW-NPTS2)/2.
  !
  !            WORK  - Real
  !                    Vector of dimension at least LENW
  !                    on return
  !                    WORK(1), ..., WORK(LAST) contain the left
  !                     end points of the subintervals in the
  !                     partition of (A,B),
  !                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
  !                     the right end points,
  !                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                     the integral approximations over the subintervals,
  !                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
  !                     contain the corresponding error estimates,
  !                    WORK(LIMIT*4+1), ..., WORK(LIMIT*4+NPTS2)
  !                     contain the integration limits and the
  !                     break points sorted in an ascending sequence.
  !                    note that LIMIT = (LENIW-NPTS2)/2.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  QAGPE, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  QAGP
  !
  REAL A , Abserr , B , Epsabs , Epsrel , F , Points , Result , Work
  INTEGER Ier , Iwork , Leniw , Lenw , limit , lvl , l1 , l2 , l3 , Neval , &
    Npts2
  !
  DIMENSION Iwork(*) , Points(*) , Work(*)
  !
  EXTERNAL F
  !
  !         CHECK VALIDITY OF LIMIT AND LENW.
  !
  !***FIRST EXECUTABLE STATEMENT  QAGP
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0.0E+00
  Abserr = 0.0E+00
  IF ( Leniw>=(3*Npts2-2).AND.Lenw>=(Leniw*2-Npts2).AND.Npts2>=2 ) THEN
    !
    !         PREPARE CALL FOR QAGPE.
    !
    limit = (Leniw-Npts2)/2
    l1 = limit + 1
    l2 = limit + l1
    l3 = limit + l2
    l4 = limit + l3
    !
    CALL QAGPE(F,A,B,Npts2,Points,Epsabs,Epsrel,limit,Result,Abserr,Neval,&
      Ier,Work(1),Work(l1),Work(l2),Work(l3),Work(l4),Iwork(1),&
      Iwork(l1),Iwork(l2),Last)
    !
    !         CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  ENDIF
  IF ( Ier==6 ) lvl = 1
  IF ( Ier/=0 ) CALL XERMSG('SLATEC','QAGP','ABNORMAL RETURN',Ier,lvl)
END SUBROUTINE QAGP
