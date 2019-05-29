!** SOSEQS
SUBROUTINE SOSEQS(FNC,N,S,Rtolx,Atolx,Tolf,Iflag,Mxit,Ncjs,Nsrrc,Nsri,&
    Iprint,Fmax,C,Nc,B,P,Temp,X,Y,Fac,Is)
  !>
  !  Subsidiary to SOS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (SOSEQS-S, DSOSEQ-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     SOSEQS solves a system of N simultaneous nonlinear equations.
  !     See the comments in the interfacing routine SOS for a more
  !     detailed description of some of the items in the calling list.
  !
  !- *******************************************************************
  !
  !   -INPUT-
  !     FNC -Function subprogram which evaluates the equations
  !     N   -Number of equations
  !     S   -Solution vector of initial guesses
  !     RTOLX-Relative error tolerance on solution components
  !     ATOLX-Absolute error tolerance on solution components
  !     TOLF-Residual error tolerance
  !     MXIT-Maximum number of allowable iterations.
  !     NCJS-Maximum number of consecutive iterative steps to perform
  !          using the same triangular Jacobian matrix approximation.
  !     NSRRC-Number of consecutive iterative steps for which the
  !          limiting precision accuracy test must be satisfied
  !          before the routine exits with IFLAG=4.
  !     NSRI-Number of consecutive iterative steps for which the
  !          diverging condition test must be satisfied before
  !          the routine exits with IFLAG=7.
  !     IPRINT-Internal printing parameter.  You must set IPRINT=-1 if you
  !          want the intermediate solution iterates and a residual norm
  !          to be printed.
  !     C   -Internal work array, dimensioned at least N*(N+1)/2.
  !     NC  -Dimension of C array. NC  .GE.  N*(N+1)/2.
  !     B   -Internal work array, dimensioned N.
  !     P   -Internal work array, dimensioned N.
  !     TEMP-Internal work array, dimensioned N.
  !     X   -Internal work array, dimensioned N.
  !     Y   -Internal work array, dimensioned N.
  !     FAC -Internal work array, dimensioned N.
  !     IS  -Internal work array, dimensioned N.
  !
  !   -OUTPUT-
  !     S   -Solution vector
  !     IFLAG-Status indicator flag
  !     MXIT-The actual number of iterations performed
  !     FMAX-Residual norm
  !     C   -Upper unit triangular matrix which approximates the
  !          forward triangularization of the full Jacobian matrix.
  !          stored in a vector with dimension at least N*(N+1)/2.
  !     B   -Contains the residuals (function values) divided
  !          by the corresponding components of the P vector
  !     P   -Array used to store the partial derivatives. After
  !          each iteration P(K) contains the maximal derivative
  !          occurring in the K-th reduced equation.
  !     TEMP-Array used to store the previous solution iterate.
  !     X   -Solution vector. Contains the values achieved on the
  !          last iteration loop upon exit from SOS.
  !     Y   -Array containing the solution increments.
  !     FAC -Array containing factors used in computing numerical
  !          derivatives.
  !     IS  -Records the pivotal information (column interchanges)
  !
  !- *********************************************************************
  !- ** Three machine dependent parameters appear in this subroutine.
  !
  !- ** The smallest positive magnitude, zero, is defined by the function
  !- ** routine R1MACH(1).
  !
  !- ** URO, The computer unit roundoff value, is defined by R1MACH(3) for
  !- ** machines that round or R1MACH(4) for machines that truncate.
  !- ** URO is the smallest positive number such that 1.+URO  .GT.  1.
  !
  !- ** The output tape unit number, LOUN, is defined by the function
  !- ** I1MACH(2).
  !- *********************************************************************
  !
  !***
  ! **See also:**  SOS
  !***
  ! **Routines called:**  I1MACH, R1MACH, SOSSOL

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  USE service, ONLY : R1MACH, I1MACH
  REAL, EXTERNAL :: FNC
  INTEGER :: Mxit, N, Nc, Ncjs, Nsri, Nsrrc, Iflag, Iprint, Is(N)
  REAL :: Atolx, Fmax, Rtolx, Tolf, C(Nc), B(N), Fac(N), P(N), S(N), Temp(N), &
    X(N), Y(N)
  INTEGER :: ksv, l, loun, ls, m, mit, mm, np1, ic, icr, isj, isv, it, item, &
    itry, j, jk, js, k, kd, kj, kk, km1, kn
  REAL :: csv, f, fact, fdif, fmin, fmxs, fn1, fn2, fp, h, hx, pmax, re, &
    sruro, test, uro, xnorm, yj, yn1, yn2, yn3, ynorm, yns, zero
  !* FIRST EXECUTABLE STATEMENT  SOSEQS
  uro = R1MACH(4)
  loun = I1MACH(2)
  zero = R1MACH(1)
  re = MAX(Rtolx,uro)
  sruro = SQRT(uro)
  !
  Iflag = 0
  np1 = N + 1
  icr = 0
  ic = 0
  itry = Ncjs
  yn1 = 0.
  yn2 = 0.
  yn3 = 0.
  yns = 0.
  mit = 0
  fn1 = 0.
  fn2 = 0.
  fmxs = 0.
  !
  !     INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND
  !     SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE.
  !
  DO k = 1, N
    Is(k) = k
    X(k) = S(k)
    Temp(k) = X(k)
  END DO
  !
  !
  !    *****************************************
  !    **** BEGIN PRINCIPAL ITERATION LOOP  ****
  !    *****************************************
  !
  DO m = 1, Mxit
    !
    DO k = 1, N
      Fac(k) = sruro
    END DO
    DO
      !
      kn = 1
      Fmax = 0.
      !
      !
      !    ******** BEGIN SUBITERATION LOOP DEFINING THE LINEARIZATION OF EACH
      !    ******** EQUATION WHICH RESULTS IN THE CONSTRUCTION OF AN UPPER
      !    ******** TRIANGULAR MATRIX APPROXIMATING THE FORWARD
      !    ******** TRIANGULARIZATION OF THE FULL JACOBIAN MATRIX
      !
      DO k = 1, N
        km1 = k - 1
        !
        !     BACK-SOLVE A TRIANGULAR LINEAR SYSTEM OBTAINING
        !     IMPROVED SOLUTION VALUES FOR K-1 OF THE VARIABLES
        !     FROM THE FIRST K-1 EQUATIONS. THESE VARIABLES ARE THEN
        !     ELIMINATED FROM THE K-TH EQUATION.
        !
        IF ( km1/=0 ) THEN
          CALL SOSSOL(k,N,km1,Y,C,B,kn)
          DO j = 1, km1
            js = Is(j)
            X(js) = Temp(js) + Y(j)
          END DO
        END IF
        !
        !
        !     EVALUATE THE K-TH EQUATION AND THE INTERMEDIATE COMPUTATION
        !     FOR THE MAX NORM OF THE RESIDUAL VECTOR.
        !
        f = FNC(X,k)
        Fmax = MAX(Fmax,ABS(f))
        !
        !     IF WE WISH TO PERFORM SEVERAL ITERATIONS USING A FIXED
        !     FACTORIZATION OF AN APPROXIMATE JACOBIAN,WE NEED ONLY
        !     UPDATE THE CONSTANT VECTOR.
        !
        IF ( itry>=Ncjs ) THEN
          !
          !
          it = 0
          !
          !     COMPUTE PARTIAL DERIVATIVES THAT ARE REQUIRED IN THE LINEARIZATION
          !     OF THE K-TH REDUCED EQUATION
          !
          DO j = k, N
            item = Is(j)
            hx = X(item)
            h = Fac(item)*hx
            IF ( ABS(h)<=zero ) h = Fac(item)
            X(item) = hx + h
            IF ( km1/=0 ) THEN
              Y(j) = h
              CALL SOSSOL(k,N,j,Y,C,B,kn)
              DO l = 1, km1
                ls = Is(l)
                X(ls) = Temp(ls) + Y(l)
              END DO
            END IF
            fp = FNC(X,k)
            X(item) = hx
            fdif = fp - f
            IF ( ABS(fdif)<=uro*ABS(f) ) THEN
              fdif = 0.
              it = it + 1
            END IF
            P(j) = fdif/h
          END DO
          !
          IF ( it>(N-k) ) THEN
            !
            !     ALL COMPUTED PARTIAL DERIVATIVES OF THE K-TH EQUATION
            !     ARE EFFECTIVELY ZERO.TRY LARGER PERTURBATIONS OF THE
            !     INDEPENDENT VARIABLES.
            !
            DO j = k, N
              isj = Is(j)
              fact = 100.*Fac(isj)
              IF ( fact>1.E+10 ) GOTO 300
              Fac(isj) = fact
            END DO
            GOTO 50
            !
          ELSEIF ( k/=N ) THEN
            !
            !     ACHIEVE A PIVOTING EFFECT BY CHOOSING THE MAXIMAL DERIVATIVE
            !     ELEMENT
            !
            pmax = 0.
            DO j = k, N
              test = ABS(P(j))
              IF ( test>pmax ) THEN
                pmax = test
                isv = j
              END IF
            END DO
            IF ( pmax==0. ) GOTO 300
            !
            !     SET UP THE COEFFICIENTS FOR THE K-TH ROW OF THE TRIANGULAR
            !     LINEAR SYSTEM AND SAVE THE PARTIAL DERIVATIVE OF
            !     LARGEST MAGNITUDE
            !
            pmax = P(isv)
            kk = kn
            DO j = k, N
              IF ( j/=isv ) C(kk) = -P(j)/pmax
              kk = kk + 1
            END DO
            P(k) = pmax
            !
            !
            IF ( isv/=k ) THEN
              !
              !     INTERCHANGE THE TWO COLUMNS OF C DETERMINED BY THE
              !     PIVOTAL STRATEGY
              !
              ksv = Is(k)
              Is(k) = Is(isv)
              Is(isv) = ksv
              !
              kd = isv - k
              kj = k
              DO j = 1, k
                csv = C(kj)
                jk = kj + kd
                C(kj) = C(jk)
                C(jk) = csv
                kj = kj + N - j
              END DO
            END IF
          END IF
        END IF
        !
        kn = kn + np1 - k
        !
        !     STORE THE COMPONENTS FOR THE CONSTANT VECTOR
        !
        B(k) = -f/P(k)
        !
      END DO
      !
      !    ********
      !    ******** END OF LOOP CREATING THE TRIANGULAR LINEARIZATION MATRIX
      !    ********
      !
      !
      !     SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW SOLUTION
      !     APPROXIMATION AND OBTAIN THE SOLUTION INCREMENT NORM.
      !
      kn = kn - 1
      Y(N) = B(N)
      IF ( N>1 ) CALL SOSSOL(N,N,N,Y,C,B,kn)
      xnorm = 0.
      ynorm = 0.
      DO j = 1, N
        yj = Y(j)
        ynorm = MAX(ynorm,ABS(yj))
        js = Is(j)
        X(js) = Temp(js) + yj
        xnorm = MAX(xnorm,ABS(X(js)))
      END DO
      !
      !
      !     PRINT INTERMEDIATE SOLUTION ITERATES AND RESIDUAL NORM IF DESIRED
      !
      IF ( Iprint==(-1) ) THEN
        mm = m - 1
        WRITE (loun,99001) Fmax, mm, (X(j),j=1,N)
        99001 FORMAT ('0RESIDUAL NORM =',E9.2,/1X,'SOLUTION ITERATE',' (',I3,')',&
          /(1X,5E26.14))
      END IF
      !
      !     TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE AND/OR ABSOLUTE ERROR
      !     COMPARISON ON SUCCESSIVE APPROXIMATIONS OF EACH SOLUTION VARIABLE)
      !
      DO j = 1, N
        js = Is(j)
        IF ( ABS(Y(j))>re*ABS(X(js))+Atolx ) GOTO 100
      END DO
      IF ( Fmax<=fmxs ) Iflag = 1
      EXIT
      50 CONTINUE
    END DO
    !
    !     TEST FOR CONVERGENCE TO A SOLUTION BASED ON RESIDUALS
    !
    100  IF ( Fmax<=Tolf ) Iflag = Iflag + 2
    IF ( Iflag>0 ) GOTO 200
    !
    !
    IF ( m>1 ) THEN
      !
      !     SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM.
      !
      IF ( Fmax<fmin ) THEN
        mit = m + 1
        yn1 = ynorm
        yn2 = yns
        fn1 = fmxs
        fmin = Fmax
        DO j = 1, N
          S(j) = X(j)
        END DO
        ic = 0
      END IF
      !
      !     TEST FOR LIMITING PRECISION CONVERGENCE.  VERY SLOWLY CONVERGENT
      !     PROBLEMS MAY ALSO BE DETECTED.
      !
      IF ( ynorm<=sruro*xnorm ) THEN
        IF ( (Fmax>=0.2*fmxs).AND.(Fmax<=5.*fmxs) ) THEN
          IF ( (ynorm>=0.2*yns).AND.(ynorm<=5.*yns) ) THEN
            icr = icr + 1
            IF ( icr<Nsrrc ) THEN
              ic = 0
              GOTO 150
            ELSE
              Iflag = 4
              Fmax = fmin
              GOTO 400
            END IF
          END IF
        END IF
      END IF
      icr = 0
      !
      !     TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME.
      !
      IF ( (ynorm<=2.*yns).AND.(Fmax<=2.*fmxs) ) THEN
        ic = 0
      ELSE
        ic = ic + 1
        IF ( ic>=Nsri ) THEN
          Iflag = 7
          GOTO 200
        END IF
      END IF
    ELSE
      fmin = Fmax
    END IF
    !
    !     CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD JACOBIAN
    !     FACTORIZATION
    !
    150  itry = itry - 1
    IF ( itry==0 ) THEN
      itry = Ncjs
    ELSEIF ( 20.*ynorm>xnorm ) THEN
      itry = Ncjs
    ELSEIF ( ynorm>2.*yns ) THEN
      itry = Ncjs
    ELSEIF ( Fmax>=2.*fmxs ) THEN
      itry = Ncjs
    END IF
    !
    !     SAVE THE CURRENT SOLUTION APPROXIMATION AND THE RESIDUAL AND
    !     SOLUTION INCREMENT NORMS FOR USE IN THE NEXT ITERATION.
    !
    DO j = 1, N
      Temp(j) = X(j)
    END DO
    IF ( m==mit ) THEN
      fn2 = Fmax
      yn3 = ynorm
    END IF
    fmxs = Fmax
    yns = ynorm
    !
    !
  END DO
  !
  !    *****************************************
  !    **** END OF PRINCIPAL ITERATION LOOP ****
  !    *****************************************
  !
  !
  !     TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED.
  m = Mxit
  Iflag = 5
  IF ( yn1>10.0*yn2.OR.yn3>10.0*yn1 ) Iflag = 6
  IF ( fn1>5.0*fmin.OR.fn2>5.0*fmin ) Iflag = 6
  IF ( Fmax>5.0*fmin ) Iflag = 6
  !
  !
  200 CONTINUE
  DO j = 1, N
    S(j) = X(j)
  END DO
  GOTO 400
  !
  !
  !     A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR.
  300 CONTINUE
  IFlag = 8
  DO j = 1, N
    S(j) = Temp(j)
  END DO
  !
  !
  400  Mxit = m
END SUBROUTINE SOSEQS
