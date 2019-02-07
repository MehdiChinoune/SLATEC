!*==DQAGS.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQAGS
SUBROUTINE DQAGS(F,A,B,Epsabs,Epsrel,Result,Abserr,Neval,Ier,Limit,Lenw,&
    Last,Iwork,Work)
  IMPLICIT NONE
  !*--DQAGS6
  !***BEGIN PROLOGUE  DQAGS
  !***PURPOSE  The routine calculates an approximation result to a given
  !            Definite integral  I = Integral of F over (A,B),
  !            Hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***LIBRARY   SLATEC (QUADPACK)
  !***CATEGORY  H2A1A1
  !***TYPE      DOUBLE PRECISION (QAGS-S, DQAGS-D)
  !***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
  !             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
  !             QUADPACK, QUADRATURE
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
  !        Double precision version
  !
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Double precision
  !                     Function subprogram defining the integrand
  !                     Function F(X). The actual name for F needs to be
  !                     Declared E X T E R N A L in the driver program.
  !
  !            A      - Double precision
  !                     Lower limit of integration
  !
  !            B      - Double precision
  !                     Upper limit of integration
  !
  !            EPSABS - Double precision
  !                     Absolute accuracy requested
  !            EPSREL - Double precision
  !                     Relative accuracy requested
  !                     If  EPSABS.LE.0
  !                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     The routine will end with IER = 6.
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
  !                             The estimates for integral and error are
  !                             less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more sub-
  !                             divisions by increasing the value of LIMIT
  !                             (and taking the according dimension
  !                             adjustments into account. However, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties. If
  !                             the position of a local difficulty can be
  !                             determined (E.G. SINGULARITY,
  !                             DISCONTINUITY WITHIN THE INTERVAL) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling the
  !                             integrator on the subranges. If possible,
  !                             an appropriate special-purpose integrator
  !                             should be used, which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 The occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                             The error may be under-estimated.
  !                         = 3 Extremely bad integrand behaviour
  !                             occurs at some points of the integration
  !                             interval.
  !                         = 4 The algorithm does not converge.
  !                             Roundoff error is detected in the
  !                             Extrapolation table. It is presumed that
  !                             the requested tolerance cannot be
  !                             achieved, and that the returned result is
  !                             the best which can be obtained.
  !                         = 5 The integral is probably divergent, or
  !                             slowly convergent. It must be noted that
  !                             divergence can occur with any other value
  !                             of IER.
  !                         = 6 The input is invalid, because
  !                             (EPSABS.LE.0 AND
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28)
  !                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
  !                             RESULT, ABSERR, NEVAL, LAST are set to
  !                             zero.  Except when LIMIT or LENW is
  !                             invalid, IWORK(1), WORK(LIMIT*2+1) and
  !                             WORK(LIMIT*3+1) are set to zero, WORK(1)
  !                             is set to A and WORK(LIMIT+1) TO B.
  !
  !         DIMENSIONING PARAMETERS
  !            LIMIT - Integer
  !                    DIMENSIONING PARAMETER FOR IWORK
  !                    LIMIT determines the maximum number of subintervals
  !                    in the partition of the given integration interval
  !                    (A,B), LIMIT.GE.1.
  !                    IF LIMIT.LT.1, the routine will end with IER = 6.
  !
  !            LENW  - Integer
  !                    DIMENSIONING PARAMETER FOR WORK
  !                    LENW must be at least LIMIT*4.
  !                    If LENW.LT.LIMIT*4, the routine will end
  !                    with IER = 6.
  !
  !            LAST  - Integer
  !                    On return, LAST equals the number of subintervals
  !                    produced in the subdivision process, determines the
  !                    number of significant elements actually in the WORK
  !                    Arrays.
  !
  !         WORK ARRAYS
  !            IWORK - Integer
  !                    Vector of dimension at least LIMIT, the first K
  !                    elements of which contain pointers
  !                    to the error estimates over the subintervals
  !                    such that WORK(LIMIT*3+IWORK(1)),... ,
  !                    WORK(LIMIT*3+IWORK(K)) form a decreasing
  !                    sequence, with K = LAST IF LAST.LE.(LIMIT/2+2),
  !                    and K = LIMIT+1-LAST otherwise
  !
  !            WORK  - Double precision
  !                    Vector of dimension at least LENW
  !                    on return
  !                    WORK(1), ..., WORK(LAST) contain the left
  !                     end-points of the subintervals in the
  !                     partition of (A,B),
  !                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
  !                     the right end-points,
  !                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                     the integral approximations over the subintervals,
  !                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
  !                     contain the error estimates.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DQAGSE, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DQAGS
  !
  !
  DOUBLE PRECISION A , Abserr , B , Epsabs , Epsrel , F , Result , Work
  INTEGER Ier , Iwork , Last , Lenw , Limit , lvl , l1 , l2 , l3 , Neval
  !
  DIMENSION Iwork(*) , Work(*)
  !
  EXTERNAL F
  !
  !         CHECK VALIDITY OF LIMIT AND LENW.
  !
  !***FIRST EXECUTABLE STATEMENT  DQAGS
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0.0D+00
  Abserr = 0.0D+00
  IF ( Limit>=1.AND.Lenw>=Limit*4 ) THEN
    !
    !         PREPARE CALL FOR DQAGSE.
    !
    l1 = Limit + 1
    l2 = Limit + l1
    l3 = Limit + l2
    !
    CALL DQAGSE(F,A,B,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,Work(1),&
      Work(l1),Work(l2),Work(l3),Iwork,Last)
    !
    !         CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  ENDIF
  IF ( Ier==6 ) lvl = 1
  IF ( Ier/=0 ) CALL XERMSG('SLATEC','DQAGS','ABNORMAL RETURN',Ier,lvl)
END SUBROUTINE DQAGS
