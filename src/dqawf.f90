!*==DQAWF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DQAWF
      SUBROUTINE DQAWF(F,A,Omega,Integr,Epsabs,Result,Abserr,Neval,Ier,Limlst,
     &                 Lst,Leniw,Maxp1,Lenw,Iwork,Work)
      IMPLICIT NONE
!*--DQAWF6
!***BEGIN PROLOGUE  DQAWF
!***PURPOSE  The routine calculates an approximation result to a given
!            Fourier integral I=Integral of F(X)*W(X) over (A,INFINITY)
!            where W(X) = COS(OMEGA*X) or W(X) = SIN(OMEGA*X).
!            Hopefully satisfying following claim for accuracy
!            ABS(I-RESULT).LE.EPSABS.
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A3A1
!***TYPE      DOUBLE PRECISION (QAWF-S, DQAWF-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, CONVERGENCE ACCELERATION,
!             FOURIER INTEGRALS, INTEGRATION BETWEEN ZEROS, QUADPACK,
!             QUADRATURE, SPECIAL-PURPOSE INTEGRAL
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Computation of Fourier integrals
!        Standard fortran subroutine
!        Double precision version
!
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
!            OMEGA  - Double precision
!                     Parameter in the integrand WEIGHT function
!
!            INTEGR - Integer
!                     Indicates which of the WEIGHT functions is used
!                     INTEGR = 1      W(X) = COS(OMEGA*X)
!                     INTEGR = 2      W(X) = SIN(OMEGA*X)
!                     IF INTEGR.NE.1.AND.INTEGR.NE.2, the routine
!                     will end with IER = 6.
!
!            EPSABS - Double precision
!                     Absolute accuracy requested, EPSABS.GT.0.
!                     If EPSABS.LE.0, the routine will end with IER = 6.
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
!                     IER.GT.0 Abnormal termination of the routine.
!                             The estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                    If OMEGA.NE.0
!                     IER = 1 Maximum number of cycles allowed
!                             has been achieved, i.e. of subintervals
!                             (A+(K-1)C,A+KC) where
!                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
!                             FOR K = 1, 2, ..., LST.
!                             One can allow more cycles by increasing
!                             the value of LIMLST (and taking the
!                             according dimension adjustments into
!                             account). Examine the array IWORK which
!                             contains the error flags on the cycles, in
!                             order to look for eventual local
!                             integration difficulties.
!                             If the position of a local difficulty
!                             can be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 4 The extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions over
!                             the cycles, does not converge to within
!                             the requested accuracy.
!                             As in the case of IER = 1, it is advised
!                             to examine the array IWORK which contains
!                             the error flags on the cycles.
!                         = 6 The input is invalid because
!                             (INTEGR.NE.1 AND INTEGR.NE.2) or
!                              EPSABS.LE.0 or LIMLST.LT.1 or
!                              LENIW.LT.(LIMLST+2) or MAXP1.LT.1 or
!                              LENW.LT.(LENIW*2+MAXP1*25).
!                              RESULT, ABSERR, NEVAL, LST are set to
!                              zero.
!                         = 7 Bad integrand behaviour occurs within
!                             one or more of the cycles. Location and
!                             type of the difficulty involved can be
!                             determined from the first LST elements of
!                             vector IWORK.  Here LST is the number of
!                             cycles actually needed (see below).
!                             IWORK(K) = 1 The maximum number of
!                                          subdivisions (=(LENIW-LIMLST)
!                                          /2) has been achieved on the
!                                          K th cycle.
!                                      = 2 Occurrence of roundoff error
!                                          is detected and prevents the
!                                          tolerance imposed on the K th
!                                          cycle, from being achieved
!                                          on this cycle.
!                                      = 3 Extremely bad integrand
!                                          behaviour occurs at some
!                                          points of the K th cycle.
!                                      = 4 The integration procedure
!                                          over the K th cycle does
!                                          not converge (to within the
!                                          required accuracy) due to
!                                          roundoff in the extrapolation
!                                          procedure invoked on this
!                                          cycle. It is assumed that the
!                                          result on this interval is
!                                          the best which can be
!                                          obtained.
!                                      = 5 The integral over the K th
!                                          cycle is probably divergent
!                                          or slowly convergent. It must
!                                          be noted that divergence can
!                                          occur with any other value of
!                                          IWORK(K).
!                    If OMEGA = 0 and INTEGR = 1,
!                    The integral is calculated by means of DQAGIE,
!                    and IER = IWORK(1) (with meaning as described
!                    for IWORK(K),K = 1).
!
!         DIMENSIONING PARAMETERS
!            LIMLST - Integer
!                     LIMLST gives an upper bound on the number of
!                     cycles, LIMLST.GE.3.
!                     If LIMLST.LT.3, the routine will end with IER = 6.
!
!            LST    - Integer
!                     On return, LST indicates the number of cycles
!                     actually needed for the integration.
!                     If OMEGA = 0, then LST is set to 1.
!
!            LENIW  - Integer
!                     Dimensioning parameter for IWORK. On entry,
!                     (LENIW-LIMLST)/2 equals the maximum number of
!                     subintervals allowed in the partition of each
!                     cycle, LENIW.GE.(LIMLST+2).
!                     If LENIW.LT.(LIMLST+2), the routine will end with
!                     IER = 6.
!
!            MAXP1  - Integer
!                     MAXP1 gives an upper bound on the number of
!                     Chebyshev moments which can be stored, i.e. for
!                     the intervals of lengths ABS(B-A)*2**(-L),
!                     L = 0,1, ..., MAXP1-2, MAXP1.GE.1.
!                     If MAXP1.LT.1, the routine will end with IER = 6.
!            LENW   - Integer
!                     Dimensioning parameter for WORK
!                     LENW must be at least LENIW*2+MAXP1*25.
!                     If LENW.LT.(LENIW*2+MAXP1*25), the routine will
!                     end with IER = 6.
!
!         WORK ARRAYS
!            IWORK  - Integer
!                     Vector of dimension at least LENIW
!                     On return, IWORK(K) FOR K = 1, 2, ..., LST
!                     contain the error flags on the cycles.
!
!            WORK   - Double precision
!                     Vector of dimension at least
!                     On return,
!                     WORK(1), ..., WORK(LST) contain the integral
!                      approximations over the cycles,
!                     WORK(LIMLST+1), ..., WORK(LIMLST+LST) contain
!                      the error estimates over the cycles.
!                     further elements of WORK have no specific
!                     meaning for the user.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DQAWFE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DQAWF
!
      DOUBLE PRECISION A , Abserr , Epsabs , F , Omega , Result , Work
      INTEGER Ier , Integr , Iwork , Leniw , Lenw , limit , Limlst , ll2 , lvl , 
     &        Lst , l1 , l2 , l3 , l4 , l5 , l6 , Maxp1 , Neval
!
      DIMENSION Iwork(*) , Work(*)
!
      EXTERNAL F
!
!         CHECK VALIDITY OF LIMLST, LENIW, MAXP1 AND LENW.
!
!***FIRST EXECUTABLE STATEMENT  DQAWF
      Ier = 6
      Neval = 0
      Result = 0.0D+00
      Abserr = 0.0D+00
      IF ( Limlst>=3.AND.Leniw>=(Limlst+2).AND.Maxp1>=1.AND.
     &     Lenw>=(Leniw*2+Maxp1*25) ) THEN
!
!         PREPARE CALL FOR DQAWFE
!
        limit = (Leniw-Limlst)/2
        l1 = Limlst + 1
        l2 = Limlst + l1
        l3 = limit + l2
        l4 = limit + l3
        l5 = limit + l4
        l6 = limit + l5
        ll2 = limit + l1
        CALL DQAWFE(F,A,Omega,Integr,Epsabs,Limlst,limit,Maxp1,Result,Abserr,
     &              Neval,Ier,Work(1),Work(l1),Iwork(1),Lst,Work(l2),Work(l3),
     &              Work(l4),Work(l5),Iwork(l1),Iwork(ll2),Work(l6))
!
!         CALL ERROR HANDLER IF NECESSARY
!
        lvl = 0
      ENDIF
      IF ( Ier==6 ) lvl = 1
      IF ( Ier/=0 ) CALL XERMSG('SLATEC','DQAWF','ABNORMAL RETURN',Ier,lvl)
      END SUBROUTINE DQAWF
