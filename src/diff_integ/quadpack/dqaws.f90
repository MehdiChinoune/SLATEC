!** DQAWS
SUBROUTINE DQAWS(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Result,Abserr,Neval,&
    Ier,Limit,Lenw,Last,Iwork,Work)
  !>
  !  The routine calculates an approximation result to a given
  !            definite integral I = Integral of F*W over (A,B),
  !            (where W shows a singular behaviour at the end points
  !            see parameter INTEGR).
  !            Hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A2A1
  !***
  ! **Type:**      DOUBLE PRECISION (QAWS-S, DQAWS-D)
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
  !                     Parameter in the integrand function, ALFA.GT.(-1)
  !                     If ALFA.LE.(-1), the routine will end with
  !                     IER = 6.
  !
  !            BETA   - Double precision
  !                     Parameter in the integrand function, BETA.GT.(-1)
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
  !         ON RETURN
  !            RESULT - Double precision
  !                     Approximation to the integral
  !
  !            ABSERR - Double precision
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
  !                     IER.GT.0 Abnormal termination of the routine
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
  !                             B.LE.A or ALFA.LE.(-1) or BETA.LE.(-1) or
  !                             or INTEGR.LT.1 or INTEGR.GT.4 or
  !                             (EPSABS.LE.0 and
  !                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
  !                             or LIMIT.LT.2 or LENW.LT.LIMIT*4.
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
  !                     integration interval (A,B), LIMIT.GE.2.
  !                     If LIMIT.LT.2, the routine will end with IER = 6.
  !
  !            LENW   - Integer
  !                     Dimensioning parameter for WORK
  !                     LENW must be at least LIMIT*4.
  !                     If LENW.LT.LIMIT*4, the routine will end
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
  !                     sequence with K = LAST if LAST.LE.(LIMIT/2+2),
  !                     and K = LIMIT+1-LAST otherwise
  !
  !            WORK   - Double precision
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
  ! **Routines called:**  DQAWSE, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : XERMSG
  !
  INTERFACE
    REAL(8) FUNCTION F(X)
      REAL(8) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Ier, Integr, Last, Lenw, Limit, Neval, Iwork(Limit)
  REAL(8) :: A, Abserr, Alfa, B, Beta, Epsabs, Epsrel, Result, Work(Lenw)
  INTEGER :: lvl, l1, l2, l3
  !
  !         CHECK VALIDITY OF LIMIT AND LENW.
  !
  !* FIRST EXECUTABLE STATEMENT  DQAWS
  Ier = 6
  Neval = 0
  Last = 0
  Result = 0.0D+00
  Abserr = 0.0D+00
  IF ( Limit>=2.AND.Lenw>=Limit*4 ) THEN
    !
    !         PREPARE CALL FOR DQAWSE.
    !
    l1 = Limit + 1
    l2 = Limit + l1
    l3 = Limit + l2
    !
    CALL DQAWSE(F,A,B,Alfa,Beta,Integr,Epsabs,Epsrel,Limit,Result,Abserr,&
      Neval,Ier,Work(1),Work(l1),Work(l2),Work(l3),Iwork,Last)
    !
    !         CALL ERROR HANDLER IF NECESSARY.
    !
    lvl = 0
  END IF
  IF ( Ier==6 ) lvl = 1
  IF ( Ier/=0 ) CALL XERMSG('SLATEC','DQAWS','ABNORMAL RETURN',Ier,lvl)
END SUBROUTINE DQAWS
