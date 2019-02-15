!DECK SNSQE
SUBROUTINE SNSQE(FCN,JAC,Iopt,N,X,Fvec,Tol,Nprint,Info,Wa,Lwa)
  IMPLICIT NONE
  REAL FCN
  INTEGER JAC
  !***BEGIN PROLOGUE  SNSQE
  !***PURPOSE  An easy-to-use code to find a zero of a system of N
  !            nonlinear functions in N variables by a modification of
  !            the Powell hybrid method.
  !***LIBRARY   SLATEC
  !***CATEGORY  F2A
  !***TYPE      SINGLE PRECISION (SNSQE-S, DNSQE-D)
  !***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM,
  !             POWELL HYBRID METHOD, ZEROS
  !***AUTHOR  Hiebert, K. L., (SNLA)
  !***DESCRIPTION
  !
  ! 1. Purpose.
  !
  !
  !       The purpose of SNSQE is to find a zero of a system of N non-
  !       linear functions in N variables by a modification of the Powell
  !       hybrid method.  This is done by using the more general nonlinear
  !       equation solver SNSQ.  The user must provide a subroutine which
  !       calculates the functions.  The user has the option of either to
  !       provide a subroutine which calculates the Jacobian or to let the
  !       code calculate it by a forward-difference approximation.  This
  !       code is the combination of the MINPACK codes (Argonne) HYBRD1
  !       and HYBRJ1.
  !
  !
  ! 2. Subroutine and Type Statements.
  !
  !       SUBROUTINE SNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
  !      *                  WA,LWA)
  !       INTEGER IOPT,N,NPRINT,INFO,LWA
  !       REAL TOL
  !       REAL X(N),FVEC(N),WA(LWA)
  !       EXTERNAL FCN,JAC
  !
  !
  ! 3. Parameters.
  !
  !       Parameters designated as input parameters must be specified on
  !       entry to SNSQE and are not changed on exit, while parameters
  !       designated as output parameters need not be specified on entry
  !       and are set to appropriate values on exit from SNSQE.
  !
  !       FCN is the name of the user-supplied subroutine which calculates
  !         the functions.  FCN must be declared in an EXTERNAL statement
  !         in the user calling program, and should be written as follows.
  !
  !         SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !         INTEGER N,IFLAG
  !         REAL X(N),FVEC(N)
  !         ----------
  !         Calculate the functions at X and
  !         return this vector in FVEC.
  !         ----------
  !         RETURN
  !         END
  !
  !         The value of IFLAG should not be changed by FCN unless the
  !         user wants to terminate execution of SNSQE.  In this case, set
  !         IFLAG to a negative integer.
  !
  !       JAC is the name of the user-supplied subroutine which calculates
  !         the Jacobian.  If IOPT=1, then JAC must be declared in an
  !         EXTERNAL statement in the user calling program, and should be
  !         written as follows.
  !
  !         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
  !         INTEGER N,LDFJAC,IFLAG
  !         REAL X(N),FVEC(N),FJAC(LDFJAC,N)
  !         ----------
  !         Calculate the Jacobian at X and return this
  !         matrix in FJAC.  FVEC contains the function
  !         values at X and should not be altered.
  !         ----------
  !         RETURN
  !         END
  !
  !         The value of IFLAG should not be changed by JAC unless the
  !         user wants to terminate execution of SNSQE.  In this case, set
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
  !       X is an array of length N.  On input, X must contain an initial
  !         estimate of the solution vector.  On output, X contains the
  !         final estimate of the solution vector.
  !
  !       FVEC is an output array of length N which contains the functions
  !         evaluated at the output X.
  !
  !       TOL is a non-negative input variable.  Termination occurs when
  !         the algorithm estimates that the relative error between X and
  !         the solution is at most TOL.  Section 4 contains more details
  !         about TOL.
  !
  !       NPRINT is an integer input variable that enables controlled
  !         printing of iterates if it is positive.  In this case, FCN is
  !         called with IFLAG = 0 at the beginning of the first iteration
  !         and every NPRINT iteration thereafter and immediately prior
  !         to return, with X and FVEC available for printing. Appropriate
  !         print statements must be added to FCN (see example). If NPRINT
  !         is not positive, no special calls of FCN with IFLAG = 0 are
  !         made.
  !
  !       INFO is an integer output variable.  If the user has terminated
  !         execution, INFO is set to the (negative) value of IFLAG.  See
  !         description of FCN and JAC. Otherwise, INFO is set as follows.
  !
  !         INFO = 0  improper input parameters.
  !
  !         INFO = 1  algorithm estimates that the relative error between
  !                   X and the solution is at most TOL.
  !
  !         INFO = 2  number of calls to FCN has reached or exceeded
  !                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
  !
  !         INFO = 3  TOL is too small.  No further improvement in the
  !                   approximate solution X is possible.
  !
  !         INFO = 4  iteration is not making good progress.
  !
  !         Sections 4 and 5 contain more details about INFO.
  !
  !       WA is a work array of length LWA.
  !
  !       LWA is a positive integer input variable not less than
  !         (3*N**2+13*N))/2.
  !
  !
  ! 4. Successful Completion.
  !
  !       The accuracy of SNSQE is controlled by the convergence parame-
  !       ter TOL.  This parameter is used in a test which makes a compar-
  !       ison between the approximation X and a solution XSOL.  SNSQE
  !       terminates when the test is satisfied.  If TOL is less than the
  !       machine precision (as defined by the function R1MACH(4)), then
  !       SNSQE attempts only to satisfy the test defined by the machine
  !       precision.  Further progress is not usually possible.  Unless
  !       high precision solutions are required, the recommended value
  !       for TOL is the square root of the machine precision.
  !
  !       The test assumes that the functions are reasonably well behaved,
  !       and, if the Jacobian is supplied by the user, that the functions
  !       and the Jacobian  coded consistently.  If these conditions
  !       are not satisfied, SNSQE may incorrectly indicate convergence.
  !       The coding of the Jacobian can be checked by the subroutine
  !       CHKDER.  If the Jacobian is coded correctly or IOPT=2, then
  !       the validity of the answer can be checked, for example, by
  !       rerunning SNSQE with a tighter tolerance.
  !
  !       Convergence Test.  If ENORM(Z) denotes the Euclidean norm of a
  !         vector Z, then this test attempts to guarantee that
  !
  !               ENORM(X-XSOL) .LE.  TOL*ENORM(XSOL).
  !
  !         If this condition is satisfied with TOL = 10**(-K), then the
  !         larger components of X have K significant decimal digits and
  !         INFO is set to 1.  There is a danger that the smaller compo-
  !         nents of X may have large relative errors, but the fast rate
  !         of convergence of SNSQE usually avoids this possibility.
  !
  !
  ! 5. Unsuccessful Completion.
  !
  !       Unsuccessful termination of SNSQE can be due to improper input
  !       parameters, arithmetic interrupts, an excessive number of func-
  !       tion evaluations, errors in the functions, or lack of good prog-
  !       ress.
  !
  !       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or
  !         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or
  !         LWA .LT. (3*N**2+13*N)/2.
  !
  !       Arithmetic Interrupts.  If these interrupts occur in the FCN
  !         subroutine during an early stage of the computation, they may
  !         be caused by an unacceptable choice of X by SNSQE.  In this
  !         case, it may be possible to remedy the situation by not evalu-
  !         ating the functions here, but instead setting the components
  !         of FVEC to numbers that exceed those in the initial FVEC.
  !
  !       Excessive Number of Function Evaluations.  If the number of
  !         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
  !         IOPT=2, then this indicates that the routine is converging
  !         very slowly as measured by the progress of FVEC, and INFO is
  !         set to 2.  This situation should be unusual because, as
  !         indicated below, lack of good progress is usually diagnosed
  !         earlier by SNSQE, causing termination with INFO = 4.
  !
  !       Errors in the Functions.  When IOPT=2, the choice of step length
  !         in the forward-difference approximation to the Jacobian
  !         assumes that the relative errors in the functions are of the
  !         order of the machine precision.  If this is not the case,
  !         SNSQE may fail (usually with INFO = 4).  The user should
  !         then either use SNSQ and set the step length or use IOPT=1
  !         and supply the Jacobian.
  !
  !       Lack of Good Progress.  SNSQE searches for a zero of the system
  !         by minimizing the sum of the squares of the functions.  In so
  !         doing, it can become trapped in a region where the minimum
  !         does not correspond to a zero of the system and, in this situ-
  !         ation, the iteration eventually fails to make good progress.
  !         In particular, this will happen if the system does not have a
  !         zero.  If the system has a zero, rerunning SNSQE from a dif-
  !         ferent starting point may be helpful.
  !
  !
  ! 6. Characteristics of the Algorithm.
  !
  !       SNSQE is a modification of the Powell hybrid method.  Two of
  !       its main characteristics involve the choice of the correction as
  !       a convex combination of the Newton and scaled gradient direc-
  !       tions, and the updating of the Jacobian by the rank-1 method of
  !       Broyden.  The choice of the correction guarantees (under reason-
  !       able conditions) global convergence for starting points far from
  !       the solution and a fast rate of convergence.  The Jacobian is
  !       calculated at the starting point by either the user-supplied
  !       subroutine or a forward-difference approximation, but it is not
  !       recalculated until the rank-1 method fails to produce satis-
  !       factory progress.
  !
  !       Timing.  The time required by SNSQE to solve a given problem
  !         depends on N, the behavior of the functions, the accuracy
  !         requested, and the starting point.  The number of arithmetic
  !         operations needed by SNSQE is about 11.5*(N**2) to process
  !         each evaluation of the functions (call to FCN) and 1.3*(N**3)
  !         to process each evaluation of the Jacobian (call to JAC,
  !         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
  !         the timing of SNSQE will be strongly influenced by the time
  !         spent in FCN and JAC.
  !
  !       Storage.  SNSQE requires (3*N**2 + 17*N)/2 single precision
  !         storage locations, in addition to the storage required by the
  !         program.  There are no internally declared storage arrays.
  !
  !
  ! 7. Example.
  !
  !       The problem is to determine the values of X(1), X(2), ..., X(9),
  !       which solve the system of tridiagonal equations
  !
  !       (3-2*X(1))*X(1)           -2*X(2)                   = -1
  !               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
  !                                   -X(8) + (3-2*X(9))*X(9) = -1
  !
  !       **********
  !
  !       PROGRAM TEST
  ! C
  ! C     Driver for SNSQE example.
  ! C
  !       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
  !       REAL TOL,FNORM
  !       REAL X(9),FVEC(9),WA(180)
  !       REAL ENORM,R1MACH
  !       EXTERNAL FCN
  !       DATA NWRITE /6/
  ! C
  !       IOPT = 2
  !       N = 9
  ! C
  ! C     The following starting values provide a rough solution.
  ! C
  !       DO 10 J = 1, 9
  !          X(J) = -1.E0
  !    10    CONTINUE
  !
  !       LWA = 180
  !       NPRINT = 0
  ! C
  ! C     Set TOL to the square root of the machine precision.
  ! C     Unless high precision solutions are required,
  ! C     this is the recommended setting.
  ! C
  !       TOL = SQRT(R1MACH(4))
  ! C
  !       CALL SNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
  !       FNORM = ENORM(N,FVEC)
  !       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
  !       STOP
  !  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
  !      *        5X,' EXIT PARAMETER',16X,I10 //
  !      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
  !       END
  !       SUBROUTINE FCN(N,X,FVEC,IFLAG)
  !       INTEGER N,IFLAG
  !       REAL X(N),FVEC(N)
  !       INTEGER K
  !       REAL ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
  !       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
  ! C
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
  !       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
  !
  !       EXIT PARAMETER                         1
  !
  !       FINAL APPROXIMATE SOLUTION
  !
  !       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
  !       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
  !       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
  !
  !***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
  !                 tions. In Numerical Methods for Nonlinear Algebraic
  !                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
  !                 1988.
  !***ROUTINES CALLED  SNSQ, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SNSQE
  INTEGER Iopt, N, Nprint, Info, Lwa
  REAL Tol
  REAL X(*), Fvec(*), Wa(Lwa)
  EXTERNAL FCN, JAC
  INTEGER index, j, lr, maxfev, ml, mode, mu, nfev, njev
  REAL epsfcn, factor, one, xtol, zero
  SAVE factor, one, zero
  DATA factor, one, zero/1.0E2, 1.0E0, 0.0E0/
  !***FIRST EXECUTABLE STATEMENT  SNSQE
  Info = 0
  !
  !     CHECK THE INPUT PARAMETERS FOR ERRORS.
  !
  IF ( Iopt>=1.AND.Iopt<=2.AND.N>0.AND.Tol>=zero.AND.Lwa>=(3*N**2+13*N)/2 )&
      THEN
    !
    !     CALL SNSQ.
    !
    maxfev = 100*(N+1)
    IF ( Iopt==2 ) maxfev = 2*maxfev
    xtol = Tol
    ml = N - 1
    mu = N - 1
    epsfcn = zero
    mode = 2
    DO j = 1, N
      Wa(j) = one
    ENDDO
    lr = (N*(N+1))/2
    index = 6*N + lr
    CALL SNSQ(FCN,JAC,Iopt,N,X,Fvec,Wa(index+1),N,xtol,maxfev,ml,mu,epsfcn,&
      Wa(1),mode,factor,Nprint,Info,nfev,njev,Wa(6*N+1),lr,Wa(N+1),&
      Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
    IF ( Info==5 ) Info = 4
  ENDIF
  IF ( Info==0 ) CALL XERMSG('SLATEC','SNSQE','INVALID INPUT PARAMETER.',2,&
    1)
  !
  !     LAST CARD OF SUBROUTINE SNSQE.
  !
END SUBROUTINE SNSQE
