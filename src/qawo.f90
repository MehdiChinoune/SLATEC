!*==QAWO.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QAWO
SUBROUTINE QAWO(F,A,B,Omega,Integr,Epsabs,Epsrel,Result,Abserr,Neval,Ier,&
    Leniw,Maxp1,Lenw,Last,Iwork,Work)
  IMPLICIT NONE
  !*--QAWO6
  !*** Start of declarations inserted by SPAG
  INTEGER Iwork, Last, Lenw, limit
  REAL Work
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  QAWO
  !***PURPOSE  Calculate an approximation to a given definite integral
  !             I = Integral of F(X)*W(X) over (A,B), where
  !                   W(X) = COS(OMEGA*X)
  !                or W(X) = SIN(OMEGA*X),
  !            hopefully satisfying the following claim for accuracy
  !                ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A2A1
  !***TYPE      SINGLE PRECISION (QAWO-S, DQAWO-D)
  !***KEYWORDS  AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD,
  !             EXTRAPOLATION, GLOBALLY ADAPTIVE,
  !             INTEGRAND WITH OSCILLATORY COS OR SIN FACTOR, QUADPACK,
  !             QUADRATURE, SPECIAL-PURPOSE
  !***AUTHOR  Piessens, Robert
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !           de Doncker, Elise
  !             Applied Mathematics and Programming Division
  !             K. U. Leuven
  !***DESCRIPTION
  !
  !        Computation of oscillatory integrals
  !        Standard fortran subroutine
  !        Real version
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Real
  !                     Function subprogram defining the function
  !                     F(X).  The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            A      - Real
  !                     Lower limit of integration
  !
  !            B      - Real
  !                     Upper limit of integration
  !
  !            OMEGA  - Real
  !                     Parameter in the integrand weight function
  !
  !            INTEGR - Integer
  !                     Indicates which of the weight functions is used
  !                     INTEGR = 1      W(X) = COS(OMEGA*X)
  !                     INTEGR = 2      W(X) = SIN(OMEGA*X)
  !                     If INTEGR.NE.1.AND.INTEGR.NE.2, the routine will
  !                     end with IER = 6.
  !
  !            EPSABS - Real
  !                     Absolute accuracy requested
  !            EPSREL - Real
  !                     Relative accuracy requested
  !                     If EPSABS.LE.0 and
  !                     EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     the routine will end with IER = 6.
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
  !                     Number of  integrand evaluations
  !
  !            IER    - Integer
  !                     IER = 0 Normal and reliable termination of the
  !                             routine. It is assumed that the requested
  !                             accuracy has been achieved.
  !                   - IER.GT.0 Abnormal termination of the routine.
  !                             The estimates for integral and error are
  !                             less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved(= LENIW/2). One can
  !                             allow more subdivisions by increasing the
  !                             value of LENIW (and taking the according
  !                             dimension adjustments into account).
  !                             However, if this yields no improvement it
  !                             is advised to analyze the integrand in
  !                             order to determine the integration
  !                             difficulties. If the position of a local
  !                             difficulty can be determined (e.g.
  !                             SINGULARITY, DISCONTINUITY within the
  !                             interval) one will probably gain from
  !                             splitting up the interval at this point
  !                             and calling the integrator on the
  !                             subranges. If possible, an appropriate
  !                             special-purpose integrator should be used
  !                             which is designed for handling the type of
  !                             difficulty involved.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                             The error may be under-estimated.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some interior points of the
  !                             integration interval.
  !                         = 4 The algorithm does not converge.
  !                             Roundoff error is detected in the
  !                             extrapolation table. It is presumed that
  !                             the requested tolerance cannot be achieved
  !                             due to roundoff in the extrapolation
  !                             table, and that the returned result is
  !                             the best which can be obtained.
  !                         = 5 The integral is probably divergent, or
  !                             slowly convergent. It must be noted that
  !                             divergence can occur with any other value
  !                             of IER.
  !                         = 6 The input is invalid, because
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or (INTEGR.NE.1 AND INTEGR.NE.2),
  !                             or LENIW.LT.2 OR MAXP1.LT.1 or
  !                             LENW.LT.LENIW*2+MAXP1*25.
  !                             RESULT, ABSERR, NEVAL, LAST are set to
  !                             zero. Except when LENIW, MAXP1 or LENW are
  !                             invalid, WORK(LIMIT*2+1), WORK(LIMIT*3+1),
  !                             IWORK(1), IWORK(LIMIT+1) are set to zero,
  !                             WORK(1) is set to A and WORK(LIMIT+1) to
  !                             B.
  !
  !         DIMENSIONING PARAMETERS
  !            LENIW  - Integer
  !                     Dimensioning parameter for IWORK.
  !                     LENIW/2 equals the maximum number of subintervals
  !                     allowed in the partition of the given integration
  !                     interval (A,B), LENIW.GE.2.
  !                     If LENIW.LT.2, the routine will end with IER = 6.
  !
  !            MAXP1  - Integer
  !                     Gives an upper bound on the number of Chebyshev
  !                     moments which can be stored, i.e. for the
  !                     intervals of lengths ABS(B-A)*2**(-L),
  !                     L=0,1, ..., MAXP1-2, MAXP1.GE.1
  !                     If MAXP1.LT.1, the routine will end with IER = 6.
  !
  !            LENW   - Integer
  !                     Dimensioning parameter for WORK
  !                     LENW must be at least LENIW*2+MAXP1*25.
  !                     If LENW.LT.(LENIW*2+MAXP1*25), the routine will
  !                     end with IER = 6.
  !
  !            LAST   - Integer
  !                     On return, LAST equals the number of subintervals
  !                     produced in the subdivision process, which
  !                     determines the number of significant elements
  !                     actually in the WORK ARRAYS.
  !
  !         WORK ARRAYS
  !            IWORK  - Integer
  !                     Vector of dimension at least LENIW
  !                     on return, the first K elements of which contain
  !                     pointers to the error estimates over the
  !                     subintervals, such that WORK(LIMIT*3+IWORK(1)), ..
  !                     WORK(LIMIT*3+IWORK(K)) form a decreasing
  !                     sequence, with LIMIT = LENW/2, and K = LAST
  !                     if LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
  !                     otherwise.
  !                     Furthermore, IWORK(LIMIT+1), ..., IWORK(LIMIT+
  !                     LAST) indicate the subdivision levels of the
  !                     subintervals, such that IWORK(LIMIT+I) = L means
  !                     that the subinterval numbered I is of length
  !                     ABS(B-A)*2**(1-L).
  !
  !            WORK   - Real
  !                     Vector of dimension at least LENW
  !                     On return
  !                     WORK(1), ..., WORK(LAST) contain the left
  !                      end points of the subintervals in the
  !                      partition of (A,B),
  !                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
  !                      the right end points,
  !                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                      the integral approximations over the
  !                      subintervals,
  !                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
  !                      contain the error estimates.
  !                     WORK(LIMIT*4+1), ..., WORK(LIMIT*4+MAXP1*25)
  !                      Provide space for storing the Chebyshev moments.
  !                     Note that LIMIT = LENW/2.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  QAWOE, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  QAWO
  !
  REAL A, Abserr, B, Epsabs, Epsrel, F, Omega, Result
  INTEGER Ier, Integr, Leniw, lvl, l1, l2, l3, l4, Maxp1, momcom, &
    Neval
  !
  DIMENSION Iwork(*), Work(*)
  !
  EXTERNAL F
  !
  !         CHECK VALIDITY OF LENIW, MAXP1 AND LENW.
  !
  !***FIRST EXECUTABLE STATEMENT  QAWO
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0.0E+00
  Abserr = 0.0E+00
  IF ( Leniw>=2.AND.Maxp1>=1.AND.Lenw>=(Leniw*2+Maxp1*25) ) THEN
    !
    !         PREPARE CALL FOR QAWOE
    !
    limit = Leniw/2
    l1 = limit + 1
    l2 = limit + l1
    l3 = limit + l2
    l4 = limit + l3
    CALL QAWOE(F,A,B,Omega,Integr,Epsabs,Epsrel,limit,1,Maxp1,Result,Abserr,&
      Neval,Ier,Last,Work(1),Work(l1),Work(l2),Work(l3),Iwork(1),&
      Iwork(l1),momcom,Work(l4))
    !
    !         CALL ERROR HANDLER IF NECESSARY
    !
    lvl = 0
  ENDIF
  IF ( Ier==6 ) lvl = 1
  IF ( Ier/=0 ) CALL XERMSG('SLATEC','QAWO','ABNORMAL RETURN',Ier,lvl)
END SUBROUTINE QAWO
