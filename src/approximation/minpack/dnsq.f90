!** DNSQ
SUBROUTINE DNSQ(FCN,JAC,Iopt,N,X,Fvec,Fjac,Ldfjac,Xtol,Maxfev,Ml,Mu,&
    Epsfcn,Diag,Mode,Factor,Nprint,Info,Nfev,Njev,R,Lr,Qtf,Wa1,Wa2,Wa3,Wa4)
  IMPLICIT NONE
  !>
  !***
  !  Find a zero of a system of a N nonlinear functions in N
  !            variables by a modification of the Powell hybrid method.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  F2A
  !***
  ! **Type:**      DOUBLE PRECISION (SNSQ-S, DNSQ-D)
  !***
  ! **Keywords:**  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS
  !***
  ! **Author:**  Hiebert, K. L. (SNLA)
  !***
  ! **Description:**
  !
  ! 1. Purpose.
  !
  !       The purpose of DNSQ is to find a zero of a system of N nonlinear
  !       functions in N variables by a modification of the Powell
  !       hybrid method.  The user must provide a subroutine which
  !       calculates the functions.  The user has the option of either to
  !       provide a subroutine which calculates the Jacobian or to let the
  !       code calculate it by a forward-difference approximation.
  !       This code is the combination of the MINPACK codes (Argonne)
  !       HYBRD and HYBRDJ.
  !
  ! 2. Subroutine and Type Statements.
  !
  !       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
  !      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
  !      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
  !       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
  !       DOUBLE PRECISION XTOL,EPSFCN,FACTOR
  !       DOUBLE PRECISION
  !       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
  !      *     WA1(N),WA2(N),WA3(N),WA4(N)
  !       EXTERNAL FCN,JAC
  !
  ! 3. Parameters.
  !
  !       Parameters designated as input parameters must be specified on
  !       entry to DNSQ and are not changed on exit, while parameters
  !       designated as output parameters need not be specified on entry
  !       and are set to appropriate values on exit from DNSQ.
  !
  !       FCN is the name of the user-supplied subroutine which calculates
  !         the functions.  FCN must be declared in an EXTERNAL statement
  !         in the user calling program, and should be written as follows.
  !
  !         SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !         INTEGER N,IFLAG
  !         DOUBLE PRECISION X(N),FVEC(N)
  !         ----------
  !         CALCULATE THE FUNCTIONS AT X AND
  !         RETURN THIS VECTOR IN FVEC.
  !         ----------
  !         RETURN
  !         END
  !
  !         The value of IFLAG should not be changed by FCN unless the
  !         user wants to terminate execution of DNSQ.  In this case set
  !         IFLAG to a negative integer.
  !
  !       JAC is the name of the user-supplied subroutine which calculates
  !         the Jacobian.  If IOPT=1, then JAC must be declared in an
  !         EXTERNAL statement in the user calling program, and should be
  !         written as follows.
  !
  !         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
  !         INTEGER N,LDFJAC,IFLAG
  !         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
  !         ----------
  !         Calculate the Jacobian at X and return this
  !         matrix in FJAC.  FVEC contains the function
  !         values at X and should not be altered.
  !         ----------
  !         RETURN
  !         END
  !
  !         The value of IFLAG should not be changed by JAC unless the
  !         user wants to terminate execution of DNSQ.  In this case set
  !         IFLAG to a negative integer.
  !
  !         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
  !
  !       IOPT is an input variable which specifies how the Jacobian will
  !         be calculated.  If IOPT=1, then the user must supply the
  !         Jacobian through the subroutine JAC.  If IOPT=2, then the
  !         code will approximate the Jacobian by forward-differencing.
  !
  !       N is a positive integer input variable set to the number of
  !         functions and variables.
  !
  !       X is an array of length N.  On input X must contain an initial
  !         estimate of the solution vector.  On output X contains the
  !         final estimate of the solution vector.
  !
  !       FVEC is an output array of length N which contains the functions
  !         evaluated at the output X.
  !
  !       FJAC is an output N by N array which contains the orthogonal
  !         matrix Q produced by the QR factorization of the final
  !         approximate Jacobian.
  !
  !       LDFJAC is a positive integer input variable not less than N
  !         which specifies the leading dimension of the array FJAC.
  !
  !       XTOL is a nonnegative input variable.  Termination occurs when
  !         the relative error between two consecutive iterates is at most
  !         XTOL.  Therefore, XTOL measures the relative error desired in
  !         the approximate solution.  Section 4 contains more details
  !         about XTOL.
  !
  !       MAXFEV is a positive integer input variable.  Termination occurs
  !         when the number of calls to FCN is at least MAXFEV by the end
  !         of an iteration.
  !
  !       ML is a nonnegative integer input variable which specifies the
  !         number of subdiagonals within the band of the Jacobian matrix.
  !         If the Jacobian is not banded or IOPT=1, set ML to at
  !         least N - 1.
  !
  !       MU is a nonnegative integer input variable which specifies the
  !         number of superdiagonals within the band of the Jacobian
  !         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
  !         least N - 1.
  !
  !       EPSFCN is an input variable used in determining a suitable step
  !         for the forward-difference approximation.  This approximation
  !         assumes that the relative errors in the functions are of the
  !         order of EPSFCN.  If EPSFCN is less than the machine
  !         precision, it is assumed that the relative errors in the
  !         functions are of the order of the machine precision.  If
  !         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
  !         argument).
  !
  !       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
  !         internally set.  If MODE = 2, DIAG must contain positive
  !         entries that serve as implicit (multiplicative) scale factors
  !         for the variables.
  !
  !       MODE is an integer input variable.  If MODE = 1, the variables
  !         will be scaled internally.  If MODE = 2, the scaling is
  !         specified by the input DIAG.  Other values of MODE are
  !         equivalent to MODE = 1.
  !
  !       FACTOR is a positive input variable used in determining the
  !         initial step bound.  This bound is set to the product of
  !         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
  !         FACTOR itself.  In most cases FACTOR should lie in the
  !         interval (.1,100.).  100. is a generally recommended value.
  !
  !       NPRINT is an integer input variable that enables controlled
  !         printing of iterates if it is positive.  In this case, FCN is
  !         called with IFLAG = 0 at the beginning of the first iteration
  !         and every NPRINT iterations thereafter and immediately prior
  !         to return, with X and FVEC available for printing. appropriate
  !         print statements must be added to FCN(see example).  If NPRINT
  !         is not positive, no special calls of FCN with IFLAG = 0 are
  !         made.
  !
  !       INFO is an integer output variable.  If the user has terminated
  !         execution, INFO is set to the (negative) value of IFLAG.  See
  !         description of FCN and JAC. Otherwise, INFO is set as follows.
  !
  !         INFO = 0  Improper input parameters.
  !
  !         INFO = 1  Relative error between two consecutive iterates is
  !                   at most XTOL.
  !
  !         INFO = 2  Number of calls to FCN has reached or exceeded
  !                   MAXFEV.
  !
  !         INFO = 3  XTOL is too small.  No further improvement in the
  !                   approximate solution X is possible.
  !
  !         INFO = 4  Iteration is not making good progress, as measured
  !                   by the improvement from the last five Jacobian
  !                   evaluations.
  !
  !         INFO = 5  Iteration is not making good progress, as measured
  !                   by the improvement from the last ten iterations.
  !
  !         Sections 4 and 5 contain more details about INFO.
  !
  !       NFEV is an integer output variable set to the number of calls to
  !         FCN.
  !
  !       NJEV is an integer output variable set to the number of calls to
  !         JAC. (If IOPT=2, then NJEV is set to zero.)
  !
  !       R is an output array of length LR which contains the upper
  !         triangular matrix produced by the QR factorization of the
  !         final approximate Jacobian, stored rowwise.
  !
  !       LR is a positive integer input variable not less than
  !         (N*(N+1))/2.
  !
  !       QTF is an output array of length N which contains the vector
  !         (Q transpose)*FVEC.
  !
  !       WA1, WA2, WA3, and WA4 are work arrays of length N.
  !
  !
  ! 4. Successful completion.
  !
  !       The accuracy of DNSQ is controlled by the convergence parameter
  !       XTOL.  This parameter is used in a test which makes a comparison
  !       between the approximation X and a solution XSOL.  DNSQ
  !       terminates when the test is satisfied.  If the convergence
  !       parameter is less than the machine precision (as defined by the
  !       function D1MACH(4)), then DNSQ only attempts to satisfy the test
  !       defined by the machine precision.  Further progress is not
  !       usually possible.
  !
  !       The test assumes that the functions are reasonably well behaved,
  !       and, if the Jacobian is supplied by the user, that the functions
  !       and the Jacobian are coded consistently.  If these conditions
  !       are not satisfied, then DNSQ may incorrectly indicate
  !       convergence.  The coding of the Jacobian can be checked by the
  !       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
  !       then the validity of the answer can be checked, for example, by
  !       rerunning DNSQ with a tighter tolerance.
  !
  !       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
  !         vector Z and D is the diagonal matrix whose entries are
  !         defined by the array DIAG, then this test attempts to
  !         guarantee that
  !
  !               DENORM(D*(X-XSOL)) .LE. XTOL*DENORM(D*XSOL).
  !
  !         If this condition is satisfied with XTOL = 10**(-K), then the
  !         larger components of D*X have K significant decimal digits and
  !         INFO is set to 1.  There is a danger that the smaller
  !         components of D*X may have large relative errors, but the fast
  !         rate of convergence of DNSQ usually avoids this possibility.
  !         Unless high precision solutions are required, the recommended
  !         value for XTOL is the square root of the machine precision.
  !
  !
  ! 5. Unsuccessful Completion.
  !
  !       Unsuccessful termination of DNSQ can be due to improper input
  !       parameters, arithmetic interrupts, an excessive number of
  !       function evaluations, or lack of good progress.
  !
  !       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
  !         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or
  !         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0,
  !         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2.
  !
  !       Arithmetic Interrupts.  If these interrupts occur in the FCN
  !         subroutine during an early stage of the computation, they may
  !         be caused by an unacceptable choice of X by DNSQ.  In this
  !         case, it may be possible to remedy the situation by rerunning
  !         DNSQ with a smaller value of FACTOR.
  !
  !       Excessive Number of Function Evaluations.  A reasonable value
  !         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
  !         If the number of calls to FCN reaches MAXFEV, then this
  !         indicates that the routine is converging very slowly as
  !         measured by the progress of FVEC, and INFO is set to 2. This
  !         situation should be unusual because, as indicated below, lack
  !         of good progress is usually diagnosed earlier by DNSQ,
  !         causing termination with info = 4 or INFO = 5.
  !
  !       Lack of Good Progress.  DNSQ searches for a zero of the system
  !         by minimizing the sum of the squares of the functions.  In so
  !         doing, it can become trapped in a region where the minimum
  !         does not correspond to a zero of the system and, in this
  !         situation, the iteration eventually fails to make good
  !         progress.  In particular, this will happen if the system does
  !         not have a zero.  If the system has a zero, rerunning DNSQ
  !         from a different starting point may be helpful.
  !
  !
  ! 6. Characteristics of The Algorithm.
  !
  !       DNSQ is a modification of the Powell Hybrid method.  Two of its
  !       main characteristics involve the choice of the correction as a
  !       convex combination of the Newton and scaled gradient directions,
  !       and the updating of the Jacobian by the rank-1 method of
  !       Broyden.  The choice of the correction guarantees (under
  !       reasonable conditions) global convergence for starting points
  !       far from the solution and a fast rate of convergence.  The
  !       Jacobian is calculated at the starting point by either the
  !       user-supplied subroutine or a forward-difference approximation,
  !       but it is not recalculated until the rank-1 method fails to
  !       produce satisfactory progress.
  !
  !       Timing.  The time required by DNSQ to solve a given problem
  !         depends on N, the behavior of the functions, the accuracy
  !         requested, and the starting point.  The number of arithmetic
  !         operations needed by DNSQ is about 11.5*(N**2) to process
  !         each evaluation of the functions (call to FCN) and 1.3*(N**3)
  !         to process each evaluation of the Jacobian (call to JAC,
  !         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
  !         the timing of DNSQ will be strongly influenced by the time
  !         spent in FCN and JAC.
  !
  !       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
  !         storage locations, in addition to the storage required by the
  !         program.  There are no internally declared storage arrays.
  !
  !- Long Description:
  !
  ! 7. Example.
  !
  !       The problem is to determine the values of X(1), X(2), ..., X(9),
  !       which solve the system of tridiagonal equations
  !
  !       (3-2*X(1))*X(1)           -2*X(2)                   = -1
  !               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
  !                                   -X(8) + (3-2*X(9))*X(9) = -1
  ! C     **********
  !
  !       PROGRAM TEST
  ! C
  ! C     Driver for DNSQ example.
  ! C
  !       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
  !      *        NWRITE
  !       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM
  !       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
  !      *     WA1(9),WA2(9),WA3(9),WA4(9)
  !       DOUBLE PRECISION DENORM,D1MACH
  !       EXTERNAL FCN
  !       DATA NWRITE /6/
  ! C
  !       IOPT = 2
  !       N = 9
  ! C
  ! C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
  ! C
  !       DO 10 J = 1, 9
  !          X(J) = -1.E0
  !    10    CONTINUE
  ! C
  !       LDFJAC = 9
  !       LR = 45
  ! C
  ! C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
  ! C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
  ! C     THIS IS THE RECOMMENDED SETTING.
  ! C
  !       XTOL = SQRT(D1MACH(4))
  ! C
  !       MAXFEV = 2000
  !       ML = 1
  !       MU = 1
  !       EPSFCN = 0.E0
  !       MODE = 2
  !       DO 20 J = 1, 9
  !          DIAG(J) = 1.E0
  !    20    CONTINUE
  !       FACTOR = 1.E2
  !       NPRINT = 0
  ! C
  !       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
  !      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
  !      *           R,LR,QTF,WA1,WA2,WA3,WA4)
  !       FNORM = DENORM(N,FVEC)
  !       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
  !       STOP
  !  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
  !      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
  !      *        5X,' EXIT PARAMETER',16X,I10 //
  !      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
  !       END
  !       SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !       INTEGER N,IFLAG
  !       DOUBLE PRECISION X(N),FVEC(N)
  !       INTEGER K
  !       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
  !       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
  ! C
  !       IF (IFLAG .NE. 0) GO TO 5
  ! C
  ! C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       DO 10 K = 1, N
  !          TEMP = (THREE - TWO*X(K))*X(K)
  !          TEMP1 = ZERO
  !          IF (K .NE. 1) TEMP1 = X(K-1)
  !          TEMP2 = ZERO
  !          IF (K .NE. N) TEMP2 = X(K+1)
  !          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
  !    10    CONTINUE
  !       RETURN
  !       END
  !
  !       Results obtained with different compilers or machines
  !       may be slightly different.
  !
  !       Final L2 norm of the residuals  0.1192636E-07
  !
  !       Number of function evaluations        14
  !
  !       Exit parameter                         1
  !
  !       Final approximate solution
  !
  !       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
  !       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
  !       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
  !
  !***
  ! **References:**  M. J. D. Powell, A hybrid method for nonlinear equa-
  !                 tions. In Numerical Methods for Nonlinear Algebraic
  !                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
  !                 1988.
  !***
  ! **Routines called:**  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1,
  !                    DQFORM, DQRFAC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  REAL(8) :: D1MACH, DENORM
  INTEGER i, iflag, Info, Iopt, iter, iwa(1), j, jm1, l, Ldfjac, &
    Lr, Maxfev, Ml, Mode, Mu, N, ncfail, ncsuc, Nfev, Njev, &
    Nprint, nslow1, nslow2
  REAL(8) :: actred, delta, Diag(*), Epsfcn, epsmch, Factor, Fjac(Ldfjac,*), &
    fnorm, fnorm1, Fvec(*), pnorm, prered, Qtf(*), R(*), ratio, summ, temp, &
    Wa1(*), Wa2(*), Wa3(*), Wa4(*), X(*), xnorm, Xtol
  EXTERNAL :: FCN
  LOGICAL jeval, sing
  REAL(8), PARAMETER :: one = 1.0D0, p1 = 1.0D-1, p5 = 5.0D-1, p001 = 1.0D-3, &
    p0001 = 1.0D-4, zero = 0.0D0
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 320
  !* FIRST EXECUTABLE STATEMENT  DNSQ
  epsmch = D1MACH(4)
  !
  Info = 0
  iflag = 0
  Nfev = 0
  Njev = 0
  !
  !        CHECK THE INPUT PARAMETERS FOR ERRORS.
  !
  !     ...EXIT
  IF ( Iopt<1.OR.Iopt>2.OR.N<=0.OR.Xtol<zero.OR.Maxfev<=0.OR.Ml<0.OR.&
    Mu<0.OR.Factor<=zero.OR.Ldfjac<N.OR.Lr<(N*(N+1))/2 ) GOTO 300
  IF ( Mode==2 ) THEN
    DO j = 1, N
      !     .........EXIT
      IF ( Diag(j)<=zero ) GOTO 300
    END DO
  END IF
  !
  !        EVALUATE THE FUNCTION AT THE STARTING POINT
  !        AND CALCULATE ITS NORM.
  !
  iflag = 1
  CALL FCN(N,X,Fvec,iflag)
  Nfev = 1
  !     ...EXIT
  IF ( iflag<0 ) GOTO 300
  fnorm = DENORM(N,Fvec)
  !
  !        INITIALIZE ITERATION COUNTER AND MONITORS.
  !
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
  !
  !        BEGINNING OF THE OUTER LOOP.
  !
  !           BEGIN BLOCK PERMITTING ...EXITS TO 90
  100  jeval = .TRUE.
  !
  !              CALCULATE THE JACOBIAN MATRIX.
  !
  IF ( Iopt==2 ) THEN
    !
    !                 CODE APPROXIMATES THE JACOBIAN
    !
    iflag = 2
    CALL DFDJC1(FCN,N,X,Fvec,Fjac,Ldfjac,iflag,Ml,Mu,Epsfcn,Wa1,Wa2)
    Nfev = Nfev + MIN(Ml+Mu+1,N)
  ELSE
    !
    !                 USER SUPPLIES JACOBIAN
    !
    CALL JAC(N,X,Fvec,Fjac,Ldfjac,iflag)
    Njev = Njev + 1
  END IF
  !
  !     .........EXIT
  IF ( iflag<0 ) GOTO 300
  !
  !              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
  !
  CALL DQRFAC(N,N,Fjac,Ldfjac,.FALSE.,iwa,1,Wa1,Wa2,Wa3)
  !
  !              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
  !              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
  !
  !           ...EXIT
  IF ( iter==1 ) THEN
    IF ( Mode/=2 ) THEN
      DO j = 1, N
        Diag(j) = Wa2(j)
        IF ( Wa2(j)==zero ) Diag(j) = one
      END DO
    END IF
    !
    !              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
    !              X AND INITIALIZE THE STEP BOUND DELTA.
    !
    DO j = 1, N
      Wa3(j) = Diag(j)*X(j)
    END DO
    xnorm = DENORM(N,Wa3)
    delta = Factor*xnorm
    IF ( delta==zero ) delta = Factor
  END IF
  !
  !           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
  !
  DO i = 1, N
    Qtf(i) = Fvec(i)
  END DO
  DO j = 1, N
    IF ( Fjac(j,j)/=zero ) THEN
      summ = zero
      DO i = j, N
        summ = summ + Fjac(i,j)*Qtf(i)
      END DO
      temp = -summ/Fjac(j,j)
      DO i = j, N
        Qtf(i) = Qtf(i) + Fjac(i,j)*temp
      END DO
    END IF
  END DO
  !
  !           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
  !
  sing = .FALSE.
  DO j = 1, N
    l = j
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO i = 1, jm1
        R(l) = Fjac(i,j)
        l = l + N - i
      END DO
    END IF
    R(l) = Wa1(j)
    IF ( Wa1(j)==zero ) sing = .TRUE.
  END DO
  !
  !           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
  !
  CALL DQFORM(N,N,Fjac,Ldfjac,Wa1)
  !
  !           RESCALE IF NECESSARY.
  !
  IF ( Mode/=2 ) THEN
    DO j = 1, N
      Diag(j) = MAX(Diag(j),Wa2(j))
    END DO
  END IF
  !
  !           BEGINNING OF THE INNER LOOP.
  !
  !
  !              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
  !
  200 CONTINUE
  IF ( Nprint>0 ) THEN
    iflag = 0
    IF ( MOD(iter-1,Nprint)==0 ) CALL FCN(N,X,Fvec,iflag)
    !     ............EXIT
    IF ( iflag<0 ) GOTO 300
  END IF
  !
  !              DETERMINE THE DIRECTION P.
  !
  CALL DDOGLG(N,R,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
  !
  !              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
  !
  DO j = 1, N
    Wa1(j) = -Wa1(j)
    Wa2(j) = X(j) + Wa1(j)
    Wa3(j) = Diag(j)*Wa1(j)
  END DO
  pnorm = DENORM(N,Wa3)
  !
  !              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
  !
  IF ( iter==1 ) delta = MIN(delta,pnorm)
  !
  !              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
  !
  iflag = 1
  CALL FCN(N,Wa2,Wa4,iflag)
  Nfev = Nfev + 1
  !     .........EXIT
  IF ( iflag>=0 ) THEN
    fnorm1 = DENORM(N,Wa4)
    !
    !              COMPUTE THE SCALED ACTUAL REDUCTION.
    !
    actred = -one
    IF ( fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
    !
    !              COMPUTE THE SCALED PREDICTED REDUCTION.
    !
    l = 1
    DO i = 1, N
      summ = zero
      DO j = i, N
        summ = summ + R(l)*Wa1(j)
        l = l + 1
      END DO
      Wa3(i) = Qtf(i) + summ
    END DO
    temp = DENORM(N,Wa3)
    prered = zero
    IF ( temp<fnorm ) prered = one - (temp/fnorm)**2
    !
    !              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
    !              REDUCTION.
    !
    ratio = zero
    IF ( prered>zero ) ratio = actred/prered
    !
    !              UPDATE THE STEP BOUND.
    !
    IF ( ratio>=p1 ) THEN
      ncfail = 0
      ncsuc = ncsuc + 1
      IF ( ratio>=p5.OR.ncsuc>1 ) delta = MAX(delta,pnorm/p5)
      IF ( ABS(ratio-one)<=p1 ) delta = pnorm/p5
    ELSE
      ncsuc = 0
      ncfail = ncfail + 1
      delta = p5*delta
    END IF
    !
    !              TEST FOR SUCCESSFUL ITERATION.
    !
    IF ( ratio>=p0001 ) THEN
      !
      !                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
      !
      DO j = 1, N
        X(j) = Wa2(j)
        Wa2(j) = Diag(j)*X(j)
        Fvec(j) = Wa4(j)
      END DO
      xnorm = DENORM(N,Wa2)
      fnorm = fnorm1
      iter = iter + 1
    END IF
    !
    !              DETERMINE THE PROGRESS OF THE ITERATION.
    !
    nslow1 = nslow1 + 1
    IF ( actred>=p001 ) nslow1 = 0
    IF ( jeval ) nslow2 = nslow2 + 1
    IF ( actred>=p1 ) nslow2 = 0
    !
    !              TEST FOR CONVERGENCE.
    !
    IF ( delta<=Xtol*xnorm.OR.fnorm==zero ) Info = 1
    !     .........EXIT
    IF ( Info==0 ) THEN
      !
      !              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
      !
      IF ( Nfev>=Maxfev ) Info = 2
      IF ( p1*MAX(p1*delta,pnorm)<=epsmch*xnorm ) Info = 3
      IF ( nslow2==5 ) Info = 4
      IF ( nslow1==10 ) Info = 5
      !     .........EXIT
      IF ( Info==0 ) THEN
        !
        !              CRITERION FOR RECALCULATING JACOBIAN
        !
        !           ...EXIT
        IF ( ncfail==2 ) GOTO 100
        !
        !              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
        !              AND UPDATE QTF IF NECESSARY.
        !
        DO j = 1, N
          summ = zero
          DO i = 1, N
            summ = summ + Fjac(i,j)*Wa4(i)
          END DO
          Wa2(j) = (summ-Wa3(j))/pnorm
          Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
          IF ( ratio>=p0001 ) Qtf(j) = summ
        END DO
        !
        !              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
        !
        CALL D1UPDT(N,N,R,Lr,Wa1,Wa2,Wa3,sing)
        CALL D1MPYQ(N,N,Fjac,Ldfjac,Wa2,Wa3)
        CALL D1MPYQ(1,N,Qtf,1,Wa2,Wa3)
        !
        !              END OF THE INNER LOOP.
        !
        jeval = .FALSE.
        !
        !           END OF THE OUTER LOOP.
        !
        GOTO 200
      END IF
    END IF
  END IF
  !
  !     TERMINATION, EITHER NORMAL OR USER IMPOSED.
  !
  300 CONTINUE
  IF ( iflag<0 ) Info = iflag
  iflag = 0
  IF ( Nprint>0 ) CALL FCN(N,X,Fvec,iflag)
  IF ( Info<0 ) CALL XERMSG('SLATEC','DNSQ',&
    'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.',1,1)
  IF ( Info==0 ) CALL XERMSG('SLATEC','DNSQ','INVALID INPUT PARAMETER.',2,1)
  IF ( Info==2 ) CALL XERMSG('SLATEC','DNSQ',&
    'TOO MANY FUNCTION EVALUATIONS.',9,1)
  IF ( Info==3 ) CALL XERMSG('SLATEC','DNSQ',&
    'XTOL TOO SMALL. NO FURTHER IMPROVEMENT POSSIBLE.',3,1)
  IF ( Info>4 ) CALL XERMSG('SLATEC','DNSQ',&
    'ITERATION NOT MAKING GOOD PROGRESS.',1,1)
  !
  !     LAST CARD OF SUBROUTINE DNSQ.
  !
END SUBROUTINE DNSQ
