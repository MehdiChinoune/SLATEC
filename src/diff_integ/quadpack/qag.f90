!** QAG
PURE SUBROUTINE QAG(F,A,B,Epsabs,Epsrel,Key,Result,Abserr,Neval,Ier,Limit,Lenw,&
    Last,Iwork,Work)
  !> The routine calculates an approximation result to a given definite integral
  !  I = integral of F over (A,B), hopefully satisfying following claim for accuracy
  !  ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A1A1
  !***
  ! **Type:**      SINGLE PRECISION (QAG-S, DQAG-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES, GENERAL-PURPOSE,
  !                GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR, QUADPACK, QUADRATURE
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
  !        Real version
  !
  !            F      - Real
  !                     Function subprogram defining the integrand
  !                     Function F(X). The actual name for F needs to be
  !                     Declared E X T E R N A L in the driver program.
  !
  !            A      - Real
  !                     Lower limit of integration
  !
  !            B      - Real
  !                     Upper limit of integration
  !
  !            EPSABS - Real
  !                     Absolute accuracy requested
  !            EPSREL - Real
  !                     Relative accuracy requested
  !                     If  EPSABS<=0
  !                     And EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     The routine will end with IER = 6.
  !
  !            KEY    - Integer
  !                     Key for choice of local integration rule
  !                     A GAUSS-KRONROD PAIR is used with
  !                       7 - 15 POINTS If KEY<2,
  !                      10 - 21 POINTS If KEY = 2,
  !                      15 - 31 POINTS If KEY = 3,
  !                      20 - 41 POINTS If KEY = 4,
  !                      25 - 51 POINTS If KEY = 5,
  !                      30 - 61 POINTS If KEY>5.
  !
  !         ON RETURN
  !            RESULT - Real
  !                     Approximation to the integral
  !
  !            ABSERR - Real
  !                     Estimate of the modulus of the absolute error,
  !                     Which should EQUAL or EXCEED ABS(I-RESULT)
  !
  !            NEVAL  - Integer
  !                     Number of integrand evaluations
  !
  !            IER    - Integer
  !                     IER = 0 Normal and reliable termination of the
  !                             routine. It is assumed that the requested
  !                             accuracy has been achieved.
  !                     IER>0 Abnormal termination of the routine
  !                             The estimates for RESULT and ERROR are
  !                             Less reliable. It is assumed that the
  !                             requested accuracy has not been achieved.
  !                      ERROR MESSAGES
  !                     IER = 1 Maximum number of subdivisions allowed
  !                             has been achieved. One can allow more
  !                             subdivisions by increasing the value of
  !                             LIMIT (and taking the according dimension
  !                             adjustments into account). HOWEVER, If
  !                             this yield no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties.
  !                             If the position of a local difficulty can
  !                             be determined (I.E. SINGULARITY,
  !                             DISCONTINUITY WITHIN THE INTERVAL) One
  !                             will probably gain from splitting up the
  !                             interval at this point and calling the
  !                             INTEGRATOR on the SUBRANGES. If possible,
  !                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
  !                             should be used which is designed for
  !                             handling the type of difficulty involved.
  !                         = 2 The occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 The input is invalid, because
  !                             (EPSABS<=0 AND
  !                              EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             OR LIMIT<1 OR LENW<LIMIT*4.
  !                             RESULT, ABSERR, NEVAL, LAST are set
  !                             to zero.
  !                             EXCEPT when LENW is invalid, IWORK(1),
  !                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
  !                             set to zero, WORK(1) is set to A and
  !                             WORK(LIMIT+1) to B.
  !
  !         DIMENSIONING PARAMETERS
  !            LIMIT - Integer
  !                    Dimensioning parameter for IWORK
  !                    Limit determines the maximum number of subintervals
  !                    in the partition of the given integration interval
  !                    (A,B), LIMIT>=1.
  !                    If LIMIT<1, the routine will end with IER = 6.
  !
  !            LENW  - Integer
  !                    Dimensioning parameter for work
  !                    LENW must be at least LIMIT*4.
  !                    IF LENW<LIMIT*4, the routine will end with
  !                    IER = 6.
  !
  !            LAST  - Integer
  !                    On return, LAST equals the number of subintervals
  !                    produced in the subdivision process, which
  !                    determines the number of significant elements
  !                    actually in the WORK ARRAYS.
  !
  !         WORK ARRAYS
  !            IWORK - Integer
  !                    Vector of dimension at least limit, the first K
  !                    elements of which contain pointers to the error
  !                    estimates over the subintervals, such that
  !                    WORK(LIMIT*3+IWORK(1)),..., WORK(LIMIT*3+IWORK(K))
  !                    form a decreasing sequence with K = LAST If
  !                    LAST<=(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
  !
  !            WORK  - Real
  !                    Vector of dimension at least LENW
  !                    on return
  !                    WORK(1), ..., WORK(LAST) contain the left end
  !                    points of the subintervals in the partition of
  !                     (A,B),
  !                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
  !                     right end points,
  !                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                     the integral approximations over the subintervals,
  !                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
  !                     the error estimates.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  QAGE, XERMSG

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
  INTEGER, INTENT(IN) :: Key, Lenw, Limit
  INTEGER, INTENT(OUT) :: Ier, Last, Neval, Iwork(Limit)
  REAL(SP), INTENT(IN) :: A, B, Epsabs, Epsrel
  REAL(SP), INTENT(OUT) :: Abserr, Result, Work(Lenw)
  !
  INTEGER :: lvl, l1, l2, l3
  !* FIRST EXECUTABLE STATEMENT  QAG
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0._SP
  Abserr = 0._SP
  IF( Limit>=1 .AND. Lenw>=Limit*4 ) THEN
    !
    !        PREPARE CALL FOR QAGE.
    !
    l1 = Limit + 1
    l2 = Limit + l1
    l3 = Limit + l2
    !
    CALL QAGE(F,A,B,Epsabs,Epsrel,Key,Limit,Result,Abserr,Neval,Ier,Work(1),&
      Work(l1),Work(l2),Work(l3),Iwork,Last)
    !
    !        CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  END IF
  !
  IF( Ier==6 ) lvl = 1
  IF( Ier/=0 ) ERROR STOP 'QAG : ABNORMAL RETURN'
  !
END SUBROUTINE QAG