!** DQAWC
SUBROUTINE DQAWC(F,A,B,C,Epsabs,Epsrel,Result,Abserr,Neval,Ier,Limit,Lenw,&
    Last,Iwork,Work)
  !>
  !  The routine calculates an approximation result to a
  !            Cauchy principal value I = INTEGRAL of F*W over (A,B)
  !            (W(X) = 1/((X-C), C.NE.A, C.NE.B), hopefully satisfying
  !            following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABE,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A1, J4
  !***
  ! **Type:**      DOUBLE PRECISION (QAWC-S, DQAWC-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE,
  !             CLENSHAW-CURTIS METHOD, GLOBALLY ADAPTIVE, QUADPACK,
  !             QUADRATURE, SPECIAL-PURPOSE
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
  !        Computation of a Cauchy principal value
  !        Standard fortran subroutine
  !        Double precision version
  !
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Double precision
  !                     Function subprogram defining the integrand
  !                     Function F(X). The actual name for F needs to be
  !                     declared E X T E R N A L in the driver program.
  !
  !            A      - Double precision
  !                     Under limit of integration
  !
  !            B      - Double precision
  !                     Upper limit of integration
  !
  !            C      - Parameter in the weight function, C.NE.A, C.NE.B.
  !                     If C = A or C = B, the routine will end with
  !                     IER = 6 .
  !
  !            EPSABS - Double precision
  !                     Absolute accuracy requested
  !            EPSREL - Double precision
  !                     Relative accuracy requested
  !                     If  EPSABS.LE.0
  !                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
  !                     the routine will end with IER = 6.
  !
  !         ON RETURN
  !            RESULT - Double precision
  !                     Approximation to the integral
  !
  !            ABSERR - Double precision
  !                     Estimate or the modulus of the absolute error,
  !                     Which should equal or exceed ABS(I-RESULT)
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
  !                             divisions by increasing the value of LIMIT
  !                             (and taking the according dimension
  !                             adjustments into account). However, if
  !                             this yields no improvement it is advised
  !                             to analyze the integrand in order to
  !                             determine the integration difficulties.
  !                             If the position of a local difficulty
  !                             can be determined (e.g. SINGULARITY,
  !                             DISCONTINUITY within the interval) one
  !                             will probably gain from splitting up the
  !                             interval at this point and calling
  !                             appropriate integrators on the subranges.
  !                         = 2 The occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 Extremely bad integrand behaviour occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 The input is invalid, because
  !                             C = A or C = B or
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or LIMIT.LT.1 or LENW.LT.LIMIT*4.
  !                             RESULT, ABSERR, NEVAL, LAST are set to
  !                             zero.  Except when LENW or LIMIT is
  !                             invalid, IWORK(1), WORK(LIMIT*2+1) and
  !                             WORK(LIMIT*3+1) are set to zero, WORK(1)
  !                             is set to A and WORK(LIMIT+1) to B.
  !
  !         DIMENSIONING PARAMETERS
  !            LIMIT - Integer
  !                    Dimensioning parameter for IWORK
  !                    LIMIT determines the maximum number of subintervals
  !                    in the partition of the given integration interval
  !                    (A,B), LIMIT.GE.1.
  !                    If LIMIT.LT.1, the routine will end with IER = 6.
  !
  !           LENW   - Integer
  !                    Dimensioning parameter for WORK
  !                    LENW must be at least LIMIT*4.
  !                    If LENW.LT.LIMIT*4, the routine will end with
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
  !                    Vector of dimension at least LIMIT, the first K
  !                    elements of which contain pointers
  !                    to the error estimates over the subintervals,
  !                    such that WORK(LIMIT*3+IWORK(1)), ... ,
  !                    WORK(LIMIT*3+IWORK(K)) form a decreasing
  !                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2),
  !                    and K = LIMIT+1-LAST otherwise
  !
  !            WORK  - Double precision
  !                    Vector of dimension at least LENW
  !                    On return
  !                    WORK(1), ..., WORK(LAST) contain the left
  !                     end points of the subintervals in the
  !                     partition of (A,B),
  !                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
  !                     the right end points,
  !                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
  !                     the integral approximations over the subintervals,
  !                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
  !                     contain the error estimates.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DQAWCE, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : XERMSG
  !
  REAL(8) :: A, Abserr, B, C, Epsabs, Epsrel, Result, Work(*)
  INTEGER Ier, Iwork(*), Last, Lenw, Limit, lvl, l1, l2, l3, Neval
  !
  REAL(8), EXTERNAL :: F
  !
  !         CHECK VALIDITY OF LIMIT AND LENW.
  !
  !* FIRST EXECUTABLE STATEMENT  DQAWC
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0.0D+00
  Abserr = 0.0D+00
  IF ( Limit>=1.AND.Lenw>=Limit*4 ) THEN
    !
    !         PREPARE CALL FOR DQAWCE.
    !
    l1 = Limit + 1
    l2 = Limit + l1
    l3 = Limit + l2
    CALL DQAWCE(F,A,B,C,Epsabs,Epsrel,Limit,Result,Abserr,Neval,Ier,Work(1),&
      Work(l1),Work(l2),Work(l3),Iwork,Last)
    !
    !         CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  END IF
  IF ( Ier==6 ) lvl = 1
  IF ( Ier/=0 ) CALL XERMSG('SLATEC','DQAWC','ABNORMAL RETURN',Ier,lvl)
END SUBROUTINE DQAWC
