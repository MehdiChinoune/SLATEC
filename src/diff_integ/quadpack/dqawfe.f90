!** DQAWFE
SUBROUTINE DQAWFE(F,A,Omega,Integr,Epsabs,Limlst,Limit,Maxp1,Result,&
    Abserr,Neval,Ier,Rslst,Erlst,Ierlst,Lst,Alist,Blist,&
    Rlist,Elist,Iord,Nnlog,Chebmo)
  !> The routine calculates an approximation result to a
  !            given Fourier integral
  !            I = Integral of F(X)*W(X) over (A,INFINITY)
  !            where W(X)=COS(OMEGA*X) or W(X)=SIN(OMEGA*X),
  !            hopefully satisfying following claim for accuracy
  !            ABS(I-RESULT)<=EPSABS.
  !***
  ! **Library:**   SLATEC (QUADPACK)
  !***
  ! **Category:**  H2A3A1
  !***
  ! **Type:**      DOUBLE PRECISION (QAWFE-S, DQAWFE-D)
  !***
  ! **Keywords:**  AUTOMATIC INTEGRATOR, CONVERGENCE ACCELERATION,
  !             FOURIER INTEGRALS, INTEGRATION BETWEEN ZEROS, QUADPACK,
  !             QUADRATURE, SPECIAL-PURPOSE INTEGRAL
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
  !        Computation of Fourier integrals
  !        Standard fortran subroutine
  !        Double precision version
  !
  !        PARAMETERS
  !         ON ENTRY
  !            F      - Double precision
  !                     Function subprogram defining the integrand
  !                     Function F(X). The actual name for F needs to
  !                     be declared E X T E R N A L in the driver program.
  !
  !            A      - Double precision
  !                     Lower limit of integration
  !
  !            OMEGA  - Double precision
  !                     Parameter in the WEIGHT function
  !
  !            INTEGR - Integer
  !                     Indicates which WEIGHT function is used
  !                     INTEGR = 1      W(X) = COS(OMEGA*X)
  !                     INTEGR = 2      W(X) = SIN(OMEGA*X)
  !                     If INTEGR/=1 .AND. INTEGR/=2, the routine will
  !                     end with IER = 6.
  !
  !            EPSABS - Double precision
  !                     absolute accuracy requested, EPSABS>0
  !                     If EPSABS<=0, the routine will end with IER = 6.
  !
  !            LIMLST - Integer
  !                     LIMLST gives an upper bound on the number of
  !                     cycles, LIMLST>=1.
  !                     If LIMLST<3, the routine will end with IER = 6.
  !
  !            LIMIT  - Integer
  !                     Gives an upper bound on the number of subintervals
  !                     allowed in the partition of each cycle, LIMIT>=1
  !                     each cycle, LIMIT>=1.
  !
  !            MAXP1  - Integer
  !                     Gives an upper bound on the number of
  !                     Chebyshev moments which can be stored, I.E.
  !                     for the intervals of lengths ABS(B-A)*2**(-L),
  !                     L=0,1, ..., MAXP1-2, MAXP1>=1
  !
  !         ON RETURN
  !            RESULT - Double precision
  !                     Approximation to the integral X
  !
  !            ABSERR - Double precision
  !                     Estimate of the modulus of the absolute error,
  !                     which should equal or exceed ABS(I-RESULT)
  !
  !            NEVAL  - Integer
  !                     Number of integrand evaluations
  !
  !            IER    - IER = 0 Normal and reliable termination of
  !                             the routine. It is assumed that the
  !                             requested accuracy has been achieved.
  !                     IER>0 Abnormal termination of the routine. The
  !                             estimates for integral and error are less
  !                             reliable. It is assumed that the requested
  !                             accuracy has not been achieved.
  !            ERROR MESSAGES
  !                    If OMEGA/=0
  !                     IER = 1 Maximum number of  cycles  allowed
  !                             Has been achieved., i.e. of subintervals
  !                             (A+(K-1)C,A+KC) where
  !                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
  !                             for K = 1, 2, ..., LST.
  !                             One can allow more cycles by increasing
  !                             the value of LIMLST (and taking the
  !                             according dimension adjustments into
  !                             account).
  !                             Examine the array IWORK which contains
  !                             the error flags on the cycles, in order to
  !                             look for eventual local integration
  !                             difficulties. If the position of a local
  !                             difficulty can be determined (e.g.
  !                             SINGULARITY, DISCONTINUITY within the
  !                             interval) one will probably gain from
  !                             splitting up the interval at this point
  !                             and calling appropriate integrators on
  !                             the subranges.
  !                         = 4 The extrapolation table constructed for
  !                             convergence acceleration of the series
  !                             formed by the integral contributions over
  !                             the cycles, does not converge to within
  !                             the requested accuracy. As in the case of
  !                             IER = 1, it is advised to examine the
  !                             array IWORK which contains the error
  !                             flags on the cycles.
  !                         = 6 The input is invalid because
  !                             (INTEGR/=1 AND INTEGR/=2) or
  !                              EPSABS<=0 or LIMLST<3.
  !                              RESULT, ABSERR, NEVAL, LST are set
  !                              to zero.
  !                         = 7 Bad integrand behaviour occurs within one
  !                             or more of the cycles. Location and type
  !                             of the difficulty involved can be
  !                             determined from the vector IERLST. Here
  !                             LST is the number of cycles actually
  !                             needed (see below).
  !                             IERLST(K) = 1 The maximum number of
  !                                           subdivisions (= LIMIT) has
  !                                           been achieved on the K th
  !                                           cycle.
  !                                       = 2 Occurrence of roundoff error
  !                                           is detected and prevents the
  !                                           tolerance imposed on the
  !                                           K th cycle, from being
  !                                           achieved.
  !                                       = 3 Extremely bad integrand
  !                                           behaviour occurs at some
  !                                           points of the K th cycle.
  !                                       = 4 The integration procedure
  !                                           over the K th cycle does
  !                                           not converge (to within the
  !                                           required accuracy) due to
  !                                           roundoff in the
  !                                           extrapolation procedure
  !                                           invoked on this cycle. It
  !                                           is assumed that the result
  !                                           on this interval is the
  !                                           best which can be obtained.
  !                                       = 5 The integral over the K th
  !                                           cycle is probably divergent
  !                                           or slowly convergent. It
  !                                           must be noted that
  !                                           divergence can occur with
  !                                           any other value of
  !                                           IERLST(K).
  !                    If OMEGA = 0 and INTEGR = 1,
  !                    The integral is calculated by means of DQAGIE
  !                    and IER = IERLST(1) (with meaning as described
  !                    for IERLST(K), K = 1).
  !
  !            RSLST  - Double precision
  !                     Vector of dimension at least LIMLST
  !                     RSLST(K) contains the integral contribution
  !                     over the interval (A+(K-1)C,A+KC) where
  !                     C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
  !                     K = 1, 2, ..., LST.
  !                     Note that, if OMEGA = 0, RSLST(1) contains
  !                     the value of the integral over (A,INFINITY).
  !
  !            ERLST  - Double precision
  !                     Vector of dimension at least LIMLST
  !                     ERLST(K) contains the error estimate corresponding
  !                     with RSLST(K).
  !
  !            IERLST - Integer
  !                     Vector of dimension at least LIMLST
  !                     IERLST(K) contains the error flag corresponding
  !                     with RSLST(K). For the meaning of the local error
  !                     flags see description of output parameter IER.
  !
  !            LST    - Integer
  !                     Number of subintervals needed for the integration
  !                     If OMEGA = 0 then LST is set to 1.
  !
  !            ALIST, BLIST, RLIST, ELIST - Double precision
  !                     vector of dimension at least LIMIT,
  !
  !            IORD, NNLOG - Integer
  !                     Vector of dimension at least LIMIT, providing
  !                     space for the quantities needed in the subdivision
  !                     process of each cycle
  !
  !            CHEBMO - Double precision
  !                     Array of dimension at least (MAXP1,25), providing
  !                     space for the Chebyshev moments needed within the
  !                     cycles
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DQAGIE, DQAWOE, DQELG

  !* REVISION HISTORY  (YYMMDD)
  !   800101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891009  Removed unreferenced variable.  (WRB)
  !   891009  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : D1MACH
  !
  INTERFACE
    REAL(DP) FUNCTION F(X)
      IMPORT DP
      REAL(DP) :: X
    END FUNCTION F
  END INTERFACE
  INTEGER :: Ier, Integr, Lst, Limit, Limlst, Maxp1, Neval
  INTEGER :: Iord(Limit), Ierlst(Limlst), Nnlog(Limit)
  REAL(DP) :: A, Abserr, Epsabs, Omega, Result
  REAL(DP) :: Alist(Limit), Blist(Limit), Chebmo(Maxp1,25), Elist(Limit), &
    Erlst(Limlst), Rlist(Limit), Rslst(Limlst)
  INTEGER :: ktmin, l, last, ll, momcom, nev, nres, numrl2
  REAL(DP) :: abseps, correc, cycle, c1, c2, dl, drl, ep, eps, epsa, errsum, &
    fact, p1, psum(52), reseps, res3la(3), uflow
  !
  !            THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF
  !            LIMEXP IN SUBROUTINE DQELG (PSUM MUST BE OF DIMENSION
  !            (LIMEXP+2) AT LEAST).
  !
  !           LIST OF MAJOR VARIABLES
  !           -----------------------
  !
  !           C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH CYCLE)
  !           CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA)
  !           PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2)
  !                       (SEE ROUTINE DQELG)
  !                       PSUM CONTAINS THE PART OF THE EPSILON TABLE
  !                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS.
  !                       EACH ELEMENT OF PSUM IS A PARTIAL SUM OF THE
  !                       SERIES WHICH SHOULD SUM TO THE VALUE OF THE
  !                       INTEGRAL.
  !           ERRSUM    - SUM OF ERROR ESTIMATES OVER THE SUBINTERVALS,
  !                       CALCULATED CUMULATIVELY
  !           EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT
  !                       SUBINTERVAL
  !           CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV
  !                       MOMENTS (SEE ALSO ROUTINE DQC25F)
  !
  REAL(DP), PARAMETER :: p = 0.9D+00
  REAL(DP), PARAMETER :: pi = 3.14159265358979323846264338327950D0
  !
  !           TEST ON VALIDITY OF PARAMETERS
  !           ------------------------------
  !
  !* FIRST EXECUTABLE STATEMENT  DQAWFE
  Result = 0.0D+00
  Abserr = 0.0D+00
  Neval = 0
  Lst = 0
  Ier = 0
  IF( (Integr/=1 .AND. Integr/=2) .OR. Epsabs<=0.0D+00 .OR. Limlst<3 ) Ier = 6
  IF( Ier/=6 ) THEN
    IF( Omega/=0.0D+00 ) THEN
      !
      !           INITIALIZATIONS
      !           ---------------
      !
      l = INT( ABS(Omega) )
      dl = 2*l + 1
      cycle = dl*pi/ABS(Omega)
      Ier = 0
      ktmin = 0
      Neval = 0
      numrl2 = 0
      nres = 0
      c1 = A
      c2 = cycle + A
      p1 = 0.1D+01 - p
      uflow = D1MACH(1)
      eps = Epsabs
      IF( Epsabs>uflow/p1 ) eps = Epsabs*p1
      ep = eps
      fact = 0.1D+01
      correc = 0.0D+00
      Abserr = 0.0D+00
      errsum = 0.0D+00
      !
      !           MAIN DO-LOOP
      !           ------------
      !
      DO Lst = 1, Limlst
        !
        !           INTEGRATE OVER CURRENT SUBINTERVAL.
        !
        epsa = eps*fact
        CALL DQAWOE(F,c1,c2,Omega,Integr,epsa,0.0D+00,Limit,Lst,Maxp1,&
          Rslst(Lst),Erlst(Lst),nev,Ierlst(Lst),last,Alist,Blist,&
          Rlist,Elist,Iord,Nnlog,momcom,Chebmo)
        Neval = Neval + nev
        fact = fact*p
        errsum = errsum + Erlst(Lst)
        drl = 0.5D+02*ABS(Rslst(Lst))
        !
        !           TEST ON ACCURACY WITH PARTIAL SUM
        !
        IF( (errsum+drl)<=Epsabs .AND. Lst>=6 ) GOTO 50
        correc = MAX(correc,Erlst(Lst))
        IF( Ierlst(Lst)/=0 ) eps = MAX(ep,correc*p1)
        IF( Ierlst(Lst)/=0 ) Ier = 7
        IF( Ier==7 .AND. (errsum+drl)<=correc*0.1D+02 .AND. Lst>5 ) GOTO 50
        numrl2 = numrl2 + 1
        IF( Lst>1 ) THEN
          psum(numrl2) = psum(ll) + Rslst(Lst)
          IF( Lst/=2 ) THEN
            !
            !           TEST ON MAXIMUM NUMBER OF SUBINTERVALS
            !
            IF( Lst==Limlst ) Ier = 1
            !
            !           PERFORM NEW EXTRAPOLATION
            !
            CALL DQELG(numrl2,psum,reseps,abseps,res3la,nres)
            !
            !           TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY ROUNDOFF
            !
            ktmin = ktmin + 1
            IF( ktmin>=15 .AND. Abserr<=0.1D-02*(errsum+drl) ) Ier = 4
            IF( abseps<=Abserr .OR. Lst==3 ) THEN
              Abserr = abseps
              Result = reseps
              ktmin = 0
              !
              !           IF IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL SUM)
              !           OR EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL
              !           APPROXIMATION
              !
              IF( (Abserr+0.1D+02*correc)<=Epsabs .OR. &
                (Abserr<=Epsabs .AND. 0.1D+02*correc>=Epsabs) ) EXIT
            END IF
            IF( Ier/=0 .AND. Ier/=7 ) EXIT
          END IF
        ELSE
          psum(1) = Rslst(1)
        END IF
        ll = numrl2
        c1 = c2
        c2 = c2 + cycle
      END DO
      !
      !         SET FINAL RESULT AND ERROR ESTIMATE
      !         -----------------------------------
      !
      Abserr = Abserr + 0.1D+02*correc
      IF( Ier==0 ) RETURN
      IF( Result==0.0D+00 .OR. psum(numrl2)==0.0D+00 ) THEN
        IF( Abserr>errsum ) GOTO 50
        IF( psum(numrl2)==0.0D+00 ) RETURN
      END IF
      IF( Abserr/ABS(Result)<=(errsum+drl)/ABS(psum(numrl2)) ) THEN
        IF( Ier>=1 .AND. Ier/=7 ) Abserr = Abserr + drl
        RETURN
      END IF
    ELSE
      !
      !           INTEGRATION BY DQAGIE IF OMEGA IS ZERO
      !           --------------------------------------
      !
      IF( Integr==1 ) CALL DQAGIE(F,A,1,Epsabs,0.0D+00,Limit,Result,Abserr,&
        Neval,Ier,Alist,Blist,Rlist,Elist,Iord,last)
      Rslst(1) = Result
      Erlst(1) = Abserr
      Ierlst(1) = Ier
      Lst = 1
      RETURN
    END IF
    50  Result = psum(numrl2)
    Abserr = errsum + drl
  END IF
  RETURN
END SUBROUTINE DQAWFE
