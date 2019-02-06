!*==DSOSEQ.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DSOSEQ
      SUBROUTINE DSOSEQ(FNC,N,S,Rtolx,Atolx,Tolf,Iflag,Mxit,Ncjs,Nsrrc,Nsri,
     &                  Iprint,Fmax,C,Nc,B,P,Temp,X,Y,Fac,Is)
      IMPLICIT NONE
!*--DSOSEQ6
!***BEGIN PROLOGUE  DSOSEQ
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DSOS
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SOSEQS-S, DSOSEQ-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     DSOSEQ solves a system of N simultaneous nonlinear equations.
!     See the comments in the interfacing routine DSOS for a more
!     detailed description of some of the items in the calling list.
!
! **********************************************************************
!   -Input-
!
!     FNC- Function subprogram which evaluates the equations
!     N  -number of equations
!     S  -Solution vector of initial guesses
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
!     IPRINT-Internal printing parameter. You must set IPRINT=-1 if you
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
!   -Output-
!     S    -Solution vector
!     IFLAG-Status indicator flag
!     MXIT-The actual number of iterations performed
!     FMAX-Residual norm
!     C   -Upper unit triangular matrix which approximates the
!          forward triangularization of the full Jacobian matrix.
!          Stored in a vector with dimension at least N*(N+1)/2.
!     B   -Contains the residuals (function values) divided
!          by the corresponding components of the P vector
!     P   -Array used to store the partial derivatives. After
!          each iteration P(K) contains the maximal derivative
!          occurring in the K-th reduced equation.
!     TEMP-Array used to store the previous solution iterate.
!     X   -Solution vector. Contains the values achieved on the
!          last iteration loop upon exit from DSOS.
!     Y   -Array containing the solution increments.
!     FAC -Array containing factors used in computing numerical
!          derivatives.
!     IS  -Records the pivotal information (column interchanges)
!
! **********************************************************************
! *** Three machine dependent parameters appear in this subroutine.
!
! *** The smallest positive magnitude, zero, is defined by the function
! *** routine D1MACH(1).
!
! *** URO, the computer unit roundoff value, is defined by D1MACH(3) for
! *** machines that round or D1MACH(4) for machines that truncate.
! *** URO is the smallest positive number such that 1.+URO  .GT.  1.
!
! *** The output tape unit number, LOUN, is defined by the function
! *** I1MACH(2).
! **********************************************************************
!
!***SEE ALSO  DSOS
!***ROUTINES CALLED  D1MACH, DSOSSL, I1MACH
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DSOSEQ
!
!
      INTEGER I1MACH
      DOUBLE PRECISION D1MACH
      INTEGER ic , icr , Iflag , Iprint , Is(*) , isj , isv , it , item , itry , 
     &        j , jk , js , k , kd , kj , kk , km1 , kn , ksv , l , loun , ls , 
     &        m , mit , mm , Mxit , N , Nc , Ncjs , np1 , Nsri , Nsrrc
      DOUBLE PRECISION Atolx , B(*) , C(*) , csv , f , Fac(*) , fact , fdif , 
     &                 Fmax , fmin , fmxs , fn1 , fn2 , FNC , fp , h , hx , 
     &                 P(*) , pmax , re , Rtolx , S(*) , sruro , Temp(*) , 
     &                 test , Tolf , uro , X(*) , xnorm , Y(*) , yj , yn1 , 
     &                 yn2 , yn3 , ynorm , yns , zero
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 430
!        BEGIN BLOCK PERMITTING ...EXITS TO 410
!           BEGIN BLOCK PERMITTING ...EXITS TO 390
!***FIRST EXECUTABLE STATEMENT  DSOSEQ
      uro = D1MACH(4)
      loun = I1MACH(2)
      zero = D1MACH(1)
      re = MAX(Rtolx,uro)
      sruro = SQRT(uro)
!
      Iflag = 0
      np1 = N + 1
      icr = 0
      ic = 0
      itry = Ncjs
      yn1 = 0.0D0
      yn2 = 0.0D0
      yn3 = 0.0D0
      yns = 0.0D0
      mit = 0
      fn1 = 0.0D0
      fn2 = 0.0D0
      fmxs = 0.0D0
!
!              INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND
!              SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE.
!
      DO k = 1 , N
        Is(k) = k
        X(k) = S(k)
        Temp(k) = X(k)
      ENDDO
!
!
!              *********************************************************
!              **** BEGIN PRINCIPAL ITERATION LOOP  ****
!              *********************************************************
!
      DO m = 1 , Mxit
!                 BEGIN BLOCK PERMITTING ...EXITS TO 350
!                    BEGIN BLOCK PERMITTING ...EXITS TO 240
!
        DO k = 1 , N
          Fac(k) = sruro
        ENDDO
        DO
!
!                          BEGIN BLOCK PERMITTING ...EXITS TO 180
          kn = 1
          Fmax = 0.0D0
!
!
!                             ******** BEGIN SUBITERATION LOOP DEFINING
!                             THE LINEARIZATION OF EACH ********
!                             EQUATION WHICH RESULTS IN THE CONSTRUCTION
!                             OF AN UPPER ******** TRIANGULAR MATRIX
!                             APPROXIMATING THE FORWARD ********
!                             TRIANGULARIZATION OF THE FULL JACOBIAN
!                             MATRIX
!
          DO k = 1 , N
!                                BEGIN BLOCK PERMITTING ...EXITS TO 160
            km1 = k - 1
!
!                                   BACK-SOLVE A TRIANGULAR LINEAR
!                                   SYSTEM OBTAINING IMPROVED SOLUTION
!                                   VALUES FOR K-1 OF THE VARIABLES FROM
!                                   THE FIRST K-1 EQUATIONS. THESE
!                                   VARIABLES ARE THEN ELIMINATED FROM
!                                   THE K-TH EQUATION.
!
            IF ( km1/=0 ) THEN
              CALL DSOSSL(k,N,km1,Y,C,B,kn)
              DO j = 1 , km1
                js = Is(j)
                X(js) = Temp(js) + Y(j)
              ENDDO
            ENDIF
!
!
!                                   EVALUATE THE K-TH EQUATION AND THE
!                                   INTERMEDIATE COMPUTATION FOR THE MAX
!                                   NORM OF THE RESIDUAL VECTOR.
!
            f = FNC(X,k)
            Fmax = MAX(Fmax,ABS(f))
!
!                                   IF WE WISH TO PERFORM SEVERAL
!                                   ITERATIONS USING A FIXED
!                                   FACTORIZATION OF AN APPROXIMATE
!                                   JACOBIAN,WE NEED ONLY UPDATE THE
!                                   CONSTANT VECTOR.
!
!                                ...EXIT
            IF ( itry>=Ncjs ) THEN
!
!
              it = 0
!
!                                   COMPUTE PARTIAL DERIVATIVES THAT ARE
!                                   REQUIRED IN THE LINEARIZATION OF THE
!                                   K-TH REDUCED EQUATION
!
              DO j = k , N
                item = Is(j)
                hx = X(item)
                h = Fac(item)*hx
                IF ( ABS(h)<=zero ) h = Fac(item)
                X(item) = hx + h
                IF ( km1/=0 ) THEN
                  Y(j) = h
                  CALL DSOSSL(k,N,j,Y,C,B,kn)
                  DO l = 1 , km1
                    ls = Is(l)
                    X(ls) = Temp(ls) + Y(l)
                  ENDDO
                ENDIF
                fp = FNC(X,k)
                X(item) = hx
                fdif = fp - f
                IF ( ABS(fdif)<=uro*ABS(f) ) THEN
                  fdif = 0.0D0
                  it = it + 1
                ENDIF
                P(j) = fdif/h
              ENDDO
!
              IF ( it>(N-k) ) THEN
!
!                                      ALL COMPUTED PARTIAL DERIVATIVES
!                                      OF THE K-TH EQUATION ARE
!                                      EFFECTIVELY ZERO.TRY LARGER
!                                      PERTURBATIONS OF THE INDEPENDENT
!                                      VARIABLES.
!
                DO j = k , N
                  isj = Is(j)
                  fact = 100.0D0*Fac(isj)
!           ..............................EXIT
                  IF ( fact>1.0D10 ) GOTO 300
                  Fac(isj) = fact
                ENDDO
!                          ............EXIT
                GOTO 50
!
!                                ...EXIT
              ELSEIF ( k/=N ) THEN
!
!                                   ACHIEVE A PIVOTING EFFECT BY
!                                   CHOOSING THE MAXIMAL DERIVATIVE
!                                   ELEMENT
!
                pmax = 0.0D0
                DO j = k , N
                  test = ABS(P(j))
                  IF ( test>pmax ) THEN
                    pmax = test
                    isv = j
                  ENDIF
                ENDDO
!           ........................EXIT
                IF ( pmax==0.0D0 ) GOTO 300
!
!                                   SET UP THE COEFFICIENTS FOR THE K-TH
!                                   ROW OF THE TRIANGULAR LINEAR SYSTEM
!                                   AND SAVE THE PARTIAL DERIVATIVE OF
!                                   LARGEST MAGNITUDE
!
                pmax = P(isv)
                kk = kn
                DO j = k , N
                  IF ( j/=isv ) C(kk) = -P(j)/pmax
                  kk = kk + 1
                ENDDO
                P(k) = pmax
!
!
!                                ...EXIT
                IF ( isv/=k ) THEN
!
!                                   INTERCHANGE THE TWO COLUMNS OF C
!                                   DETERMINED BY THE PIVOTAL STRATEGY
!
                  ksv = Is(k)
                  Is(k) = Is(isv)
                  Is(isv) = ksv
!
                  kd = isv - k
                  kj = k
                  DO j = 1 , k
                    csv = C(kj)
                    jk = kj + kd
                    C(kj) = C(jk)
                    C(jk) = csv
                    kj = kj + N - j
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
!
            kn = kn + np1 - k
!
!                                STORE THE COMPONENTS FOR THE CONSTANT
!                                VECTOR
!
            B(k) = -f/P(k)
!
!                       ......EXIT
          ENDDO
!
!                       ********
!                       ******** END OF LOOP CREATING THE TRIANGULAR
!                       LINEARIZATION MATRIX
!                       ********
!
!
!                        SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW
!                        SOLUTION APPROXIMATION AND OBTAIN THE SOLUTION
!                        INCREMENT NORM.
!
          kn = kn - 1
          Y(N) = B(N)
          IF ( N>1 ) CALL DSOSSL(N,N,N,Y,C,B,kn)
          xnorm = 0.0D0
          ynorm = 0.0D0
          DO j = 1 , N
            yj = Y(j)
            ynorm = MAX(ynorm,ABS(yj))
            js = Is(j)
            X(js) = Temp(js) + yj
            xnorm = MAX(xnorm,ABS(X(js)))
          ENDDO
!
!
!                       PRINT INTERMEDIATE SOLUTION ITERATES AND
!                       RESIDUAL NORM IF DESIRED
!
          IF ( Iprint==(-1) ) THEN
            mm = m - 1
            WRITE (loun,99001) Fmax , mm , (X(j),j=1,N)
99001       FORMAT ('0RESIDUAL NORM =',D9.2,/1X,'SOLUTION ITERATE (',I3,')',
     &              /(1X,5D26.14))
          ENDIF
!
!                       TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE
!                       AND/OR ABSOLUTE ERROR COMPARISON ON SUCCESSIVE
!                       APPROXIMATIONS OF EACH SOLUTION VARIABLE)
!
          DO j = 1 , N
            js = Is(j)
!                    ......EXIT
            IF ( ABS(Y(j))>re*ABS(X(js))+Atolx ) GOTO 100
          ENDDO
          IF ( Fmax<=fmxs ) Iflag = 1
          EXIT
 50     ENDDO
!
!                    TEST FOR CONVERGENCE TO A SOLUTION BASED ON
!                    RESIDUALS
!
 100    IF ( Fmax<=Tolf ) Iflag = Iflag + 2
!        ............EXIT
        IF ( Iflag>0 ) GOTO 200
!
!
        IF ( m>1 ) THEN
!                       BEGIN BLOCK PERMITTING ...EXITS TO 320
!
!                          SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM.
!
          IF ( Fmax<fmin ) THEN
            mit = m + 1
            yn1 = ynorm
            yn2 = yns
            fn1 = fmxs
            fmin = Fmax
            DO j = 1 , N
              S(j) = X(j)
            ENDDO
            ic = 0
          ENDIF
!
!                          TEST FOR LIMITING PRECISION CONVERGENCE. VERY
!                          SLOWLY CONVERGENT PROBLEMS MAY ALSO BE
!                          DETECTED.
!
          IF ( ynorm<=sruro*xnorm ) THEN
            IF ( Fmax>=0.2D0*fmxs.AND.Fmax<=5.0D0*fmxs ) THEN
              IF ( ynorm>=0.2D0*yns.AND.ynorm<=5.0D0*yns ) THEN
                icr = icr + 1
                IF ( icr>=Nsrrc ) THEN
                  Iflag = 4
                  Fmax = fmin
!     ........................EXIT
                  GOTO 400
                ELSE
                  ic = 0
!                       .........EXIT
                  GOTO 150
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          icr = 0
!
!                          TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME.
!
          IF ( ynorm>2.0D0*yns.OR.Fmax>2.0D0*fmxs ) THEN
            ic = ic + 1
!                       ......EXIT
            IF ( ic>=Nsri ) THEN
              Iflag = 7
!        .....................EXIT
              GOTO 200
            ENDIF
          ELSE
            ic = 0
          ENDIF
        ELSE
          fmin = Fmax
        ENDIF
!
!                    CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD
!                    JACOBIAN FACTORIZATION
!
 150    itry = itry - 1
        IF ( itry==0 ) THEN
          itry = Ncjs
        ELSEIF ( 20.0D0*ynorm>xnorm ) THEN
          itry = Ncjs
        ELSEIF ( ynorm>2.0D0*yns ) THEN
          itry = Ncjs
!                 ......EXIT
        ELSEIF ( Fmax>=2.0D0*fmxs ) THEN
          itry = Ncjs
        ENDIF
!
!                 SAVE THE CURRENT SOLUTION APPROXIMATION AND THE
!                 RESIDUAL AND SOLUTION INCREMENT NORMS FOR USE IN THE
!                 NEXT ITERATION.
!
        DO j = 1 , N
          Temp(j) = X(j)
        ENDDO
        IF ( m==mit ) THEN
          fn2 = Fmax
          yn3 = ynorm
        ENDIF
        fmxs = Fmax
        yns = ynorm
!
!
      ENDDO
!
!              *********************************************************
!              **** END OF PRINCIPAL ITERATION LOOP ****
!              *********************************************************
!
!
!               TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED.
      m = Mxit
      Iflag = 5
      IF ( yn1>10.0D0*yn2.OR.yn3>10.0D0*yn1 ) Iflag = 6
      IF ( fn1>5.0D0*fmin.OR.fn2>5.0D0*fmin ) Iflag = 6
!        ......EXIT
      IF ( Fmax>5.0D0*fmin ) Iflag = 6
!
!
 200  DO j = 1 , N
        S(j) = X(j)
      ENDDO
      GOTO 400
!
!
!           A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR.
 300  Iflag = 8
      DO j = 1 , N
        S(j) = Temp(j)
!     ......EXIT
      ENDDO
!
!
 400  Mxit = m
      END SUBROUTINE DSOSEQ
