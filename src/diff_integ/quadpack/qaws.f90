!** QAWS
PURE SUBROUTINE QAWS(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Result,Abserr,Neval,&
    Ier,Limit,Lenw,Last,Iwork,Work)
  !> The routine calculates an approximation result to a given definite integral
  !  I = Integral of F*W over (A,B),
  !  (where W shows a singular behaviour at the end points see parameter INTEGR).
  !  Hopefully satisfying following claim for accuracy
  !  ABS(I-RESULT)<=MAX(EPSABS,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A1
  !***
  ! **Type:**      SINGLE PRECISION (QAWS-S, DQAWS-D)
  !***
  ! **Keywords:**  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES,
  !             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD,
  !             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE, SPECIAL-PURPOSE
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
  !        Integration of functions having algebraico-logarithmic
  !        end point singularities
  !        Standard fortran subroutine
  !        Real version
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Real
  !                     Function subprogram defining the integrand
  !                     function F(X). The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            A      - Real
  !                     Lower limit of integration
  !
  !            B      - Real
  !                     Upper limit of integration, B>A
  !                     If B<=A, the routine will end with IER = 6.
  !
  !            ALFA   - Real
  !                     Parameter in the integrand function, ALFA>(-1)
  !                     If ALFA<=(-1), the routine will end with
  !                     IER = 6.
  !
  !            BETA   - Real
  !                     Parameter in the integrand function, BETA>(-1)
  !                     If BETA<=(-1), the routine will end with
  !                     IER = 6.
  !
  !            INTEGR - Integer
  !                     Indicates which WEIGHT function is to be used
  !                     = 1  (X-A)**ALFA*(B-X)**BETA
  !                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
  !                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
  !                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
  !                     If INTEGR<1 or INTEGR>4, the routine
  !                     will end with IER = 6.
  !
  !            EPSABS - Real
  !                     Absolute accuracy requested
  !            EPSREL - Real
  !                     Relative accuracy requested
  !                     If  EPSABS<=0
  !                     and EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     the routine will end with IER = 6.
  !
  !         ON RETURN
  !            RESULT - Real
  !                     Approximation to the integral
  !
  !            ABSERR - Real
  !                     Estimate of the modulus of the absolute error,
  !                     Which should equal or exceed ABS(I-RESULT)
  !
  !            NEVAL  - Integer
  !                     Number of integrand evaluations
  !
  !            IER    - Integer
  !                     IER = 0 Normal and reliable termination of the
  !                             routine. It is assumed that the requested
  !                             accuracy has been achieved.
  !                     IER>0 Abnormal termination of the routine
  !                             The estimates for the integral and error
  !                             are less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !            ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more
  !                             subdivisions by increasing the value of
  !                             LIMIT (and taking the according dimension
  !                             adjustments into account). However, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand, in order to
  !                             determine the integration difficulties
  !                             which prevent the requested tolerance from
  !                             being achieved. In case of a jump
  !                             discontinuity or a local singularity
  !                             of algebraico-logarithmic type at one or
  !                             more interior points of the integration
  !                             range, one should proceed by splitting up
  !                             the interval at these points and calling
  !                             the integrator on the subranges.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 The input is invalid, because
  !                             B<=A or ALFA<=(-1) or BETA<=(-1) or
  !                             or INTEGR<1 or INTEGR>4 or
  !                             (EPSABS<=0 and
  !                              EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or LIMIT<2 or LENW<LIMIT*4.
  !                             RESULT, ABSERR, NEVAL, LAST are set to
  !                             zero. Except when LENW or LIMIT is invalid
  !                             IWORK(1), WORK(LIMIT*2+1) and
  !                             WORK(LIMIT*3+1) are set to zero, WORK(1)
  !                             is set to A and WORK(LIMIT+1) to B.
  !
  !         DIMENSIONING PARAMETERS
  !            LIMIT  - Integer
  !                     Dimensioning parameter for IWORK
  !                     LIMIT determines the maximum number of
  !                     subintervals in the partition of the given
  !                     integration interval (A,B), LIMIT>=2.
  !                     If LIMIT<2, the routine will end with IER = 6.
  !
  !            LENW   - Integer
  !                     Dimensioning parameter for WORK
  !                     LENW must be at least LIMIT*4.
  !                     If LENW<LIMIT*4, the routine will end
  !                     with IER = 6.
  !
  !            LAST   - Integer
  !                     On return, LAST equals the number of
  !                     subintervals produced in the subdivision process,
  !                     which determines the significant number of
  !                     elements actually in the WORK ARRAYS.
  !
  !         WORK ARRAYS
  !            IWORK  - Integer
  !                     Vector of dimension LIMIT, the first K
  !                     elements of which contain pointers
  !                     to the error estimates over the subintervals,
  !                     such that WORK(LIMIT*3+IWORK(1)), ...,
  !                     WORK(LIMIT*3+IWORK(K)) form a decreasing
  !                     sequence with K = LAST if LAST<=(LIMIT/2+2),
  !                     and K = LIMIT+1-LAST otherwise
  !
  !            WORK   - Real
  !                     Vector of dimension LENW
  !                     On return
  !                     WORK(1), ..., WORK(LAST) contain the left
  !                      end points of the subintervals in the
  !                      partition of (A,B),
  !                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
  !                      the right end points,
  !                     WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST)
  !                      contain the integral approximations over
  !                      the subintervals,
  !                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
  !                      contain the error estimates.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  QAWSE, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  INTERFACE
    REAL(SP) PURE FUNCTION F(X)
      IMPORT SP
      REAL(SP), INTENT(IN) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER, INTENT(IN) :: Integr, Lenw, Limit
  INTEGER, INTENT(OUT) :: Ier, Last, Neval, Iwork(Limit)
  REAL(SP), INTENT(IN) :: A, Alfa, B, Beta, Epsabs, Epsrel
  REAL(SP), INTENT(OUT) :: Abserr, Result, Work(Lenw)
  !
  INTEGER :: lvl, l1, l2, l3
  !
  !         CHECK VALIDITY OF LIMIT AND LENW.
  !
  !* FIRST EXECUTABLE STATEMENT  QAWS
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0._SP
  Abserr = 0._SP
  IF( Limit>=2 .AND. Lenw>=Limit*4 ) THEN
    !
    !         PREPARE CALL FOR QAWSE.
    !
    l1 = Limit + 1
    l2 = Limit + l1
    l3 = Limit + l2
    !
    CALL QAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit,Result,Abserr,&
      Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
    !
    !         CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  END IF
  IF( Ier==6 ) lvl = 1
  IF( Ier/=0 ) ERROR STOP 'QAWS : ABNORMAL RETURN'
  !
END SUBROUTINE QAWS