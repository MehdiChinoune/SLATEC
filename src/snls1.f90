!*==SNLS1.f90  processed by SPAG 6.72Dc at 11:02 on  6 Feb 2019
!DECK SNLS1
SUBROUTINE SNLS1(FCN,Iopt,M,N,X,Fvec,Fjac,Ldfjac,Ftol,Xtol,Gtol,Maxfev,&
    Epsfcn,Diag,Mode,Factor,Nprint,Info,Nfev,Njev,Ipvt,Qtf,&
    Wa1,Wa2,Wa3,Wa4)
  IMPLICIT NONE
  !*--SNLS17
  !***BEGIN PROLOGUE  SNLS1
  !***PURPOSE  Minimize the sum of the squares of M nonlinear functions
  !            in N variables by a modification of the Levenberg-Marquardt
  !            algorithm.
  !***LIBRARY   SLATEC
  !***CATEGORY  K1B1A1, K1B1A2
  !***TYPE      SINGLE PRECISION (SNLS1-S, DNLS1-D)
  !***KEYWORDS  LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
  !             NONLINEAR LEAST SQUARES
  !***AUTHOR  Hiebert, K. L., (SNLA)
  !***DESCRIPTION
  !
  ! 1. Purpose.
  !
  !       The purpose of SNLS1 is to minimize the sum of the squares of M
  !       nonlinear functions in N variables by a modification of the
  !       Levenberg-Marquardt algorithm.  The user must provide a subrou-
  !       tine which calculates the functions.  The user has the option
  !       of how the Jacobian will be supplied.  The user can supply the
  !       full Jacobian, or the rows of the Jacobian (to avoid storing
  !       the full Jacobian), or let the code approximate the Jacobian by
  !       forward-differencing.   This code is the combination of the
  !       MINPACK codes (Argonne) LMDER, LMDIF, and LMSTR.
  !
  !
  ! 2. Subroutine and Type Statements.
  !
  !       SUBROUTINE SNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
  !      *                 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO
  !      *                 ,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
  !       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
  !       INTEGER IPVT(N)
  !       REAL FTOL,XTOL,GTOL,EPSFCN,FACTOR
  !       REAL X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),
  !      *     WA1(N),WA2(N),WA3(N),WA4(M)
  !
  !
  ! 3. Parameters.
  !
  !       Parameters designated as input parameters must be specified on
  !       entry to SNLS1 and are not changed on exit, while parameters
  !       designated as output parameters need not be specified on entry
  !       and are set to appropriate values on exit from SNLS1.
  !
  !       FCN is the name of the user-supplied subroutine which calculates
  !         the functions.  If the user wants to supply the Jacobian
  !         (IOPT=2 or 3), then FCN must be written to calculate the
  !         Jacobian, as well as the functions.  See the explanation
  !         of the IOPT argument below.
  !         If the user wants the iterates printed (NPRINT positive), then
  !         FCN must do the printing.  See the explanation of NPRINT
  !         below.  FCN must be declared in an EXTERNAL statement in the
  !         calling program and should be written as follows.
  !
  !
  !         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
  !         INTEGER IFLAG,LDFJAC,M,N
  !         REAL X(N),FVEC(M)
  !         ----------
  !         FJAC and LDFJAC may be ignored     , if IOPT=1.
  !         REAL FJAC(LDFJAC,N)                , if IOPT=2.
  !         REAL FJAC(N)                       , if IOPT=3.
  !         ----------
  !           If IFLAG=0, the values in X and FVEC are available
  !           for printing.  See the explanation of NPRINT below.
  !           IFLAG will never be zero unless NPRINT is positive.
  !           The values of X and FVEC must not be changed.
  !         RETURN
  !         ----------
  !           If IFLAG=1, calculate the functions at X and return
  !           this vector in FVEC.
  !         RETURN
  !         ----------
  !           If IFLAG=2, calculate the full Jacobian at X and return
  !           this matrix in FJAC.  Note that IFLAG will never be 2 unless
  !           IOPT=2.  FVEC contains the function values at X and must
  !           not be altered.  FJAC(I,J) must be set to the derivative
  !           of FVEC(I) with respect to X(J).
  !         RETURN
  !         ----------
  !           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
  !           and return this vector in FJAC.  Note that IFLAG will
  !           never be 3 unless IOPT=3.  FVEC contains the function
  !           values at X and must not be altered.  FJAC(J) must be
  !           set to the derivative of FVEC(LDFJAC) with respect to X(J).
  !         RETURN
  !         ----------
  !         END
  !
  !
  !         The value of IFLAG should not be changed by FCN unless the
  !         user wants to terminate execution of SNLS1.  In this case, set
  !         IFLAG to a negative integer.
  !
  !
  !       IOPT is an input variable which specifies how the Jacobian will
  !         be calculated.  If IOPT=2 or 3, then the user must supply the
  !         Jacobian, as well as the function values, through the
  !         subroutine FCN.  If IOPT=2, the user supplies the full
  !         Jacobian with one call to FCN.  If IOPT=3, the user supplies
  !         one row of the Jacobian with each call.  (In this manner,
  !         storage can be saved because the full Jacobian is not stored.)
  !         If IOPT=1, the code will approximate the Jacobian by forward
  !         differencing.
  !
  !       M is a positive integer input variable set to the number of
  !         functions.
  !
  !       N is a positive integer input variable set to the number of
  !         variables.  N must not exceed M.
  !
  !       X is an array of length N.  On input, X must contain an initial
  !         estimate of the solution vector.  On output, X contains the
  !         final estimate of the solution vector.
  !
  !       FVEC is an output array of length M which contains the functions
  !         evaluated at the output X.
  !
  !       FJAC is an output array.  For IOPT=1 and 2, FJAC is an M by N
  !         array.  For IOPT=3, FJAC is an N by N array.  The upper N by N
  !         submatrix of FJAC contains an upper triangular matrix R with
  !         diagonal elements of nonincreasing magnitude such that
  !
  !                T     T           T
  !               P *(JAC *JAC)*P = R *R,
  !
  !         where P is a permutation matrix and JAC is the final calcu-
  !         lated Jacobian.  Column J of P is column IPVT(J) (see below)
  !         of the identity matrix.  The lower part of FJAC contains
  !         information generated during the computation of R.
  !
  !       LDFJAC is a positive integer input variable which specifies
  !         the leading dimension of the array FJAC.  For IOPT=1 and 2,
  !         LDFJAC must not be less than M.  For IOPT=3, LDFJAC must not
  !         be less than N.
  !
  !       FTOL is a non-negative input variable.  Termination occurs when
  !         both the actual and predicted relative reductions in the sum
  !         of squares are at most FTOL.  Therefore, FTOL measures the
  !         relative error desired in the sum of squares.  Section 4 con-
  !         tains more details about FTOL.
  !
  !       XTOL is a non-negative input variable.  Termination occurs when
  !         the relative error between two consecutive iterates is at most
  !         XTOL.  Therefore, XTOL measures the relative error desired in
  !         the approximate solution.  Section 4 contains more details
  !         about XTOL.
  !
  !       GTOL is a non-negative input variable.  Termination occurs when
  !         the cosine of the angle between FVEC and any column of the
  !         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
  !         measures the orthogonality desired between the function vector
  !         and the columns of the Jacobian.  Section 4 contains more
  !         details about GTOL.
  !
  !       MAXFEV is a positive integer input variable.  Termination occurs
  !         when the number of calls to FCN to evaluate the functions
  !         has reached MAXFEV.
  !
  !       EPSFCN is an input variable used in determining a suitable step
  !         for the forward-difference approximation.  This approximation
  !         assumes that the relative errors in the functions are of the
  !         order of EPSFCN.  If EPSFCN is less than the machine preci-
  !         sion, it is assumed that the relative errors in the functions
  !         are of the order of the machine precision.  If IOPT=2 or 3,
  !         then EPSFCN can be ignored (treat it as a dummy argument).
  !
  !       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
  !         internally set.  If MODE = 2, DIAG must contain positive
  !         entries that serve as implicit (multiplicative) scale factors
  !         for the variables.
  !
  !       MODE is an integer input variable.  If MODE = 1, the variables
  !         will be scaled internally.  If MODE = 2, the scaling is speci-
  !         fied by the input DIAG.  Other values of MODE are equivalent
  !         to MODE = 1.
  !
  !       FACTOR is a positive input variable used in determining the ini-
  !         tial step bound.  This bound is set to the product of FACTOR
  !         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
  !         itself.  In most cases FACTOR should lie in the interval
  !         (.1,100.).  100. is a generally recommended value.
  !
  !       NPRINT is an integer input variable that enables controlled
  !         printing of iterates if it is positive.  In this case, FCN is
  !         called with IFLAG = 0 at the beginning of the first iteration
  !         and every NPRINT iterations thereafter and immediately prior
  !         to return, with X and FVEC available for printing. Appropriate
  !         print statements must be added to FCN (see example) and
  !         FVEC should not be altered.  If NPRINT is not positive, no
  !         special calls to FCN with IFLAG = 0 are made.
  !
  !       INFO is an integer output variable.  If the user has terminated
  !         execution, INFO is set to the (negative) value of IFLAG.  See
  !         description of FCN and JAC. Otherwise, INFO is set as follows.
  !
  !         INFO = 0  improper input parameters.
  !
  !         INFO = 1  both actual and predicted relative reductions in the
  !                   sum of squares are at most FTOL.
  !
  !         INFO = 2  relative error between two consecutive iterates is
  !                   at most XTOL.
  !
  !         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
  !
  !         INFO = 4  the cosine of the angle between FVEC and any column
  !                   of the Jacobian is at most GTOL in absolute value.
  !
  !         INFO = 5  number of calls to FCN for function evaluation
  !                   has reached MAXFEV.
  !
  !         INFO = 6  FTOL is too small.  No further reduction in the sum
  !                   of squares is possible.
  !
  !         INFO = 7  XTOL is too small.  No further improvement in the
  !                   approximate solution X is possible.
  !
  !         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
  !                   columns of the Jacobian to machine precision.
  !
  !         Sections 4 and 5 contain more details about INFO.
  !
  !       NFEV is an integer output variable set to the number of calls to
  !         FCN for function evaluation.
  !
  !       NJEV is an integer output variable set to the number of
  !         evaluations of the full Jacobian.  If IOPT=2, only one call to
  !         FCN is required for each evaluation of the full Jacobian.
  !         If IOPT=3, the M calls to FCN are required.
  !         If IOPT=1, then NJEV is set to zero.
  !
  !       IPVT is an integer output array of length N.  IPVT defines a
  !         permutation matrix P such that JAC*P = Q*R, where JAC is the
  !         final calculated Jacobian, Q is orthogonal (not stored), and R
  !         is upper triangular with diagonal elements of nonincreasing
  !         magnitude.  Column J of P is column IPVT(J) of the identity
  !         matrix.
  !
  !       QTF is an output array of length N which contains the first N
  !         elements of the vector (Q transpose)*FVEC.
  !
  !       WA1, WA2, and WA3 are work arrays of length N.
  !
  !       WA4 is a work array of length M.
  !
  !
  ! 4. Successful Completion.
  !
  !       The accuracy of SNLS1 is controlled by the convergence parame-
  !       ters FTOL, XTOL, and GTOL.  These parameters are used in tests
  !       which make three types of comparisons between the approximation
  !       X and a solution XSOL.  SNLS1 terminates when any of the tests
  !       is satisfied.  If any of the convergence parameters is less than
  !       the machine precision (as defined by the function R1MACH(4)),
  !       then SNLS1 only attempts to satisfy the test defined by the
  !       machine precision.  Further progress is not usually possible.
  !
  !       The tests assume that the functions are reasonably well behaved,
  !       and, if the Jacobian is supplied by the user, that the functions
  !       and the Jacobian are coded consistently.  If these conditions
  !       are not satisfied, then SNLS1 may incorrectly indicate conver-
  !       gence.  If the Jacobian is coded correctly or IOPT=1,
  !       then the validity of the answer can be checked, for example, by
  !       rerunning SNLS1 with tighter tolerances.
  !
  !       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
  !         of a vector Z, then this test attempts to guarantee that
  !
  !               ENORM(FVEC) .LE. (1+FTOL)*ENORM(FVECS),
  !
  !         where FVECS denotes the functions evaluated at XSOL.  If this
  !         condition is satisfied with FTOL = 10**(-K), then the final
  !         residual norm ENORM(FVEC) has K significant decimal digits and
  !         INFO is set to 1 (or to 3 if the second test is also satis-
  !         fied).  Unless high precision solutions are required, the
  !         recommended value for FTOL is the square root of the machine
  !         precision.
  !
  !       Second Convergence Test.  If D is the diagonal matrix whose
  !         entries are defined by the array DIAG, then this test attempts
  !         to guarantee that
  !
  !               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL).
  !
  !         If this condition is satisfied with XTOL = 10**(-K), then the
  !         larger components of D*X have K significant decimal digits and
  !         INFO is set to 2 (or to 3 if the first test is also satis-
  !         fied).  There is a danger that the smaller components of D*X
  !         may have large relative errors, but if MODE = 1, then the
  !         accuracy of the components of X is usually related to their
  !         sensitivity.  Unless high precision solutions are required,
  !         the recommended value for XTOL is the square root of the
  !         machine precision.
  !
  !       Third Convergence Test.  This test is satisfied when the cosine
  !         of the angle between FVEC and any column of the Jacobian at X
  !         is at most GTOL in absolute value.  There is no clear rela-
  !         tionship between this test and the accuracy of SNLS1, and
  !         furthermore, the test is equally well satisfied at other crit-
  !         ical points, namely maximizers and saddle points.  Therefore,
  !         termination caused by this test (INFO = 4) should be examined
  !         carefully.  The recommended value for GTOL is zero.
  !
  !
  ! 5. Unsuccessful Completion.
  !
  !       Unsuccessful termination of SNLS1 can be due to improper input
  !       parameters, arithmetic interrupts, or an excessive number of
  !       function evaluations.
  !
  !       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
  !         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or for IOPT=1 or 2
  !         LDFJAC .LT. M, or for IOPT=3 LDFJAC .LT. N, or FTOL .LT. 0.E0,
  !         or XTOL .LT. 0.E0, or GTOL .LT. 0.E0, or MAXFEV .LE. 0, or
  !         FACTOR .LE. 0.E0.
  !
  !       Arithmetic Interrupts.  If these interrupts occur in the FCN
  !         subroutine during an early stage of the computation, they may
  !         be caused by an unacceptable choice of X by SNLS1.  In this
  !         case, it may be possible to remedy the situation by rerunning
  !         SNLS1 with a smaller value of FACTOR.
  !
  !       Excessive Number of Function Evaluations.  A reasonable value
  !         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for
  !         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then
  !         this indicates that the routine is converging very slowly
  !         as measured by the progress of FVEC, and INFO is set to 5.
  !         In this case, it may be helpful to restart SNLS1 with MODE
  !         set to 1.
  !
  !
  ! 6. Characteristics of the Algorithm.
  !
  !       SNLS1 is a modification of the Levenberg-Marquardt algorithm.
  !       Two of its main characteristics involve the proper use of
  !       implicitly scaled variables (if MODE = 1) and an optimal choice
  !       for the correction.  The use of implicitly scaled variables
  !       achieves scale invariance of SNLS1 and limits the size of the
  !       correction in any direction where the functions are changing
  !       rapidly.  The optimal choice of the correction guarantees (under
  !       reasonable conditions) global convergence from starting points
  !       far from the solution and a fast rate of convergence for
  !       problems with small residuals.
  !
  !       Timing.  The time required by SNLS1 to solve a given problem
  !         depends on M and N, the behavior of the functions, the accu-
  !         racy requested, and the starting point.  The number of arith-
  !         metic operations needed by SNLS1 is about N**3 to process each
  !         evaluation of the functions (call to FCN) and to process each
  !         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one
  !         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and
  !         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
  !         can be evaluated quickly, the timing of SNLS1 will be
  !         strongly influenced by the time spent in FCN.
  !
  !       Storage.  SNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
  !         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
  !         locations and N integer storage locations, in addition to
  !         the storage required by the program.  There are no internally
  !         declared storage arrays.
  !
  ! *Long Description:
  !
  ! 7. Example.
  !
  !       The problem is to determine the values of X(1), X(2), and X(3)
  !       which provide the best fit (in the least squares sense) of
  !
  !             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
  !
  !       to the data
  !
  !             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
  !                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
  !
  !       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
  !       I-th component of FVEC is thus defined by
  !
  !             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
  !
  !       **********
  !
  !       PROGRAM TEST
  ! C
  ! C     Driver for SNLS1 example.
  ! C
  !       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,
  !      *        NWRITE
  !       INTEGER IPVT(3)
  !       REAL FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
  !       REAL X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3),
  !      *     WA1(3),WA2(3),WA3(3),WA4(15)
  !       REAL ENORM,R1MACH
  !       EXTERNAL FCN
  !       DATA NWRITE /6/
  ! C
  !       IOPT = 1
  !       M = 15
  !       N = 3
  ! C
  ! C     The following starting values provide a rough fit.
  ! C
  !       X(1) = 1.E0
  !       X(2) = 1.E0
  !       X(3) = 1.E0
  ! C
  !       LDFJAC = 15
  ! C
  ! C     Set FTOL and XTOL to the square root of the machine precision
  ! C     and GTOL to zero.  Unless high precision solutions are
  ! C     required, these are the recommended settings.
  ! C
  !       FTOL = SQRT(R1MACH(4))
  !       XTOL = SQRT(R1MACH(4))
  !       GTOL = 0.E0
  ! C
  !       MAXFEV = 400
  !       EPSFCN = 0.0
  !       MODE = 1
  !       FACTOR = 1.E2
  !       NPRINT = 0
  ! C
  !       CALL SNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
  !      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
  !      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
  !       FNORM = ENORM(M,FVEC)
  !       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
  !       STOP
  !  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
  !      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
  !      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 //
  !      *        5X,' EXIT PARAMETER',16X,I10 //
  !      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
  !       END
  !       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
  ! C     This is the form of the FCN routine if IOPT=1,
  ! C     that is, if the user does not calculate the Jacobian.
  !       INTEGER M,N,IFLAG
  !       REAL X(N),FVEC(M)
  !       INTEGER I
  !       REAL TMP1,TMP2,TMP3,TMP4
  !       REAL Y(15)
  !       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
  !      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
  !      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
  !      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
  ! C
  !       IF (IFLAG .NE. 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF (I .GT. 8) TMP3 = TMP2
  !          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
  !    10    CONTINUE
  !       RETURN
  !       END
  !
  !
  !       Results obtained with different compilers or machines
  !       may be slightly different.
  !
  !       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
  !
  !       NUMBER OF FUNCTION EVALUATIONS        25
  !
  !       NUMBER OF JACOBIAN EVALUATIONS         0
  !
  !       EXIT PARAMETER                         1
  !
  !       FINAL APPROXIMATE SOLUTION
  !
  !        0.8241058E-01  0.1133037E+01  0.2343695E+01
  !
  !
  !       For IOPT=2, FCN would be modified as follows to also
  !       calculate the full Jacobian when IFLAG=2.
  !
  !       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
  ! C
  ! C     This is the form of the FCN routine if IOPT=2,
  ! C     that is, if the user calculates the full Jacobian.
  ! C
  !       INTEGER LDFJAC,M,N,IFLAG
  !       REAL X(N),FVEC(M)
  !       REAL FJAC(LDFJAC,N)
  !       INTEGER I
  !       REAL TMP1,TMP2,TMP3,TMP4
  !       REAL Y(15)
  !       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
  !      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
  !      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
  !      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
  ! C
  !       IF (IFLAG .NE. 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       IF(IFLAG.NE.1) GO TO 20
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF (I .GT. 8) TMP3 = TMP2
  !          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
  !    10    CONTINUE
  !       RETURN
  ! C
  ! C     Below, calculate the full Jacobian.
  ! C
  !    20    CONTINUE
  ! C
  !       DO 30 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF (I .GT. 8) TMP3 = TMP2
  !          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
  !          FJAC(I,1) = -1.E0
  !          FJAC(I,2) = TMP1*TMP2/TMP4
  !          FJAC(I,3) = TMP1*TMP3/TMP4
  !    30    CONTINUE
  !       RETURN
  !       END
  !
  !
  !       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
  !         LDFJAC would be set to 3, and FCN would be written as
  !         follows to calculate a row of the Jacobian when IFLAG=3.
  !
  !       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
  ! C     This is the form of the FCN routine if IOPT=3,
  ! C     that is, if the user calculates the Jacobian row by row.
  !       INTEGER M,N,IFLAG
  !       REAL X(N),FVEC(M)
  !       REAL FJAC(N)
  !       INTEGER I
  !       REAL TMP1,TMP2,TMP3,TMP4
  !       REAL Y(15)
  !       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
  !      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
  !      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
  !      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
  ! C
  !       IF (IFLAG .NE. 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       IF( IFLAG.NE.1) GO TO 20
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF (I .GT. 8) TMP3 = TMP2
  !          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
  !    10    CONTINUE
  !       RETURN
  ! C
  ! C     Below, calculate the LDFJAC-th row of the Jacobian.
  ! C
  !    20 CONTINUE
  !
  !       I = LDFJAC
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF (I .GT. 8) TMP3 = TMP2
  !          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
  !          FJAC(1) = -1.E0
  !          FJAC(2) = TMP1*TMP2/TMP4
  !          FJAC(3) = TMP1*TMP3/TMP4
  !       RETURN
  !       END
  !
  !***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
  !                 implementation and theory.  In Numerical Analysis
  !                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
  !                 Watson, Editor), Lecture Notes in Mathematics 630,
  !                 Springer-Verlag, 1978.
  !***ROUTINES CALLED  CHKDER, ENORM, FDJAC3, LMPAR, QRFAC, R1MACH,
  !                    RWUPDT, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SNLS1
  INTEGER Iopt , M , N , Ldfjac , Maxfev , Mode , Nprint , Info , Nfev ,&
    Njev
  INTEGER ijunk , nrow , Ipvt(*)
  REAL Ftol , Xtol , Gtol , Factor , Epsfcn
  REAL X(*) , Fvec(*) , Fjac(Ldfjac,*) , Diag(*) , Qtf(*) , Wa1(*) , Wa2(*)&
    , Wa3(*) , Wa4(*)
  LOGICAL sing
  EXTERNAL FCN
  INTEGER i , iflag , iter , j , l , modech
  REAL actred , delta , dirder , epsmch , fnorm , fnorm1 , gnorm , one ,&
    par , pnorm , prered , p1 , p5 , p25 , p75 , p0001 , ratio , sum ,&
    temp , temp1 , temp2 , xnorm , zero
  REAL R1MACH , ENORM , err , chklim
  CHARACTER(8) :: xern1
  CHARACTER(16) :: xern3
  !
  SAVE chklim , one , p1 , p5 , p25 , p75 , p0001 , zero
  DATA chklim/.1E0/
  DATA one , p1 , p5 , p25 , p75 , p0001 , zero/1.0E0 , 1.0E-1 , 5.0E-1 ,&
    2.5E-1 , 7.5E-1 , 1.0E-4 , 0.0E0/
  !
  !***FIRST EXECUTABLE STATEMENT  SNLS1
  epsmch = R1MACH(4)
  !
  Info = 0
  iflag = 0
  Nfev = 0
  Njev = 0
  !
  !     CHECK THE INPUT PARAMETERS FOR ERRORS.
  !
  IF ( Iopt<1.OR.Iopt>3.OR.N<=0.OR.M<N.OR.Ldfjac<N.OR.Ftol<zero.OR.&
    Xtol<zero.OR.Gtol<zero.OR.Maxfev<=0.OR.Factor<=zero ) GOTO 200
  IF ( Iopt<3.AND.Ldfjac<M ) GOTO 200
  IF ( Mode==2 ) THEN
    DO j = 1 , N
      IF ( Diag(j)<=zero ) GOTO 200
    ENDDO
  ENDIF
  !
  !     EVALUATE THE FUNCTION AT THE STARTING POINT
  !     AND CALCULATE ITS NORM.
  !
  iflag = 1
  ijunk = 1
  CALL FCN(iflag,M,N,X,Fvec,Fjac,ijunk)
  Nfev = 1
  IF ( iflag<0 ) GOTO 200
  fnorm = ENORM(M,Fvec)
  !
  !     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
  !
  par = zero
  iter = 1
  !
  !     BEGINNING OF THE OUTER LOOP.
  !
  !
  !        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
  !
  100  IF ( Nprint>0 ) THEN
  iflag = 0
  IF ( MOD(iter-1,Nprint)==0 ) CALL FCN(iflag,M,N,X,Fvec,Fjac,ijunk)
  IF ( iflag<0 ) GOTO 200
ENDIF
!
!        CALCULATE THE JACOBIAN MATRIX.
!
IF ( Iopt==3 ) THEN
  !
  !        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
  !        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
  !        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
  !        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
  !        N COMPONENTS IN QTF.
  !
  DO j = 1 , N
    Qtf(j) = zero
    DO i = 1 , N
      Fjac(i,j) = zero
    ENDDO
  ENDDO
  DO i = 1 , M
    nrow = i
    iflag = 3
    CALL FCN(iflag,M,N,X,Fvec,Wa3,nrow)
    IF ( iflag<0 ) GOTO 200
    !
    !            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
    !
    IF ( iter<=1 ) THEN
      !
      !            GET THE INCREMENTED X-VALUES INTO WA1(*).
      !
      modech = 1
      CALL CHKDER(M,N,X,Fvec,Fjac,Ldfjac,Wa1,Wa4,modech,err)
      !
      !            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED.
      !
      IF ( i==1 ) THEN
        !
        !            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
        !
        iflag = 1
        CALL FCN(iflag,M,N,Wa1,Wa4,Fjac,nrow)
        Nfev = Nfev + 1
        IF ( iflag<0 ) GOTO 200
      ENDIF
      modech = 2
      CALL CHKDER(1,N,X,Fvec(i),Wa3,1,Wa1,Wa4(i),modech,err)
      IF ( err<chklim ) THEN
        WRITE (xern1,'(I8)') i
        WRITE (xern3,'(1PE15.6)') err
        CALL XERMSG('SLATEC','SNLS1','DERIVATIVE OF FUNCTION '//xern1//&
          ' MAY BE WRONG, ERR = '//xern3//' TOO CLOSE TO 0.',7,&
          0)
      ENDIF
    ENDIF
    !
    temp = Fvec(i)
    CALL RWUPDT(N,Fjac,Ldfjac,Wa3,Qtf,temp,Wa1,Wa2)
  ENDDO
  Njev = Njev + 1
  !
  !        IF THE JACOBIAN IS RANK DEFICIENT, CALL QRFAC TO
  !        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
  !
  sing = .FALSE.
  DO j = 1 , N
    IF ( Fjac(j,j)==zero ) sing = .TRUE.
    Ipvt(j) = j
    Wa2(j) = ENORM(j,Fjac(1,j))
  ENDDO
  IF ( sing ) THEN
    CALL QRFAC(N,N,Fjac,Ldfjac,.TRUE.,Ipvt,N,Wa1,Wa2,Wa3)
    DO j = 1 , N
      IF ( Fjac(j,j)/=zero ) THEN
        sum = zero
        DO i = j , N
          sum = sum + Fjac(i,j)*Qtf(i)
        ENDDO
        temp = -sum/Fjac(j,j)
        DO i = j , N
          Qtf(i) = Qtf(i) + Fjac(i,j)*temp
        ENDDO
      ENDIF
      Fjac(j,j) = Wa1(j)
    ENDDO
  ENDIF
ELSE
  !
  !     STORE THE FULL JACOBIAN USING M*N STORAGE
  !
  IF ( Iopt==1 ) THEN
    !
    !     THE CODE APPROXIMATES THE JACOBIAN
    !
    iflag = 1
    CALL FDJAC3(FCN,M,N,X,Fvec,Fjac,Ldfjac,iflag,Epsfcn,Wa4)
    Nfev = Nfev + N
  ELSE
    !
    !     THE USER SUPPLIES THE JACOBIAN
    !
    iflag = 2
    CALL FCN(iflag,M,N,X,Fvec,Fjac,Ldfjac)
    Njev = Njev + 1
    !
    !             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
    !
    IF ( iter<=1 ) THEN
      IF ( iflag<0 ) GOTO 200
      !
      !           GET THE INCREMENTED X-VALUES INTO WA1(*).
      !
      modech = 1
      CALL CHKDER(M,N,X,Fvec,Fjac,Ldfjac,Wa1,Wa4,modech,err)
      !
      !           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
      !
      iflag = 1
      CALL FCN(iflag,M,N,Wa1,Wa4,Fjac,Ldfjac)
      Nfev = Nfev + 1
      IF ( iflag<0 ) GOTO 200
      DO i = 1 , M
        modech = 2
        CALL CHKDER(1,N,X,Fvec(i),Fjac(i,1),Ldfjac,Wa1,Wa4(i),modech,err)
        IF ( err<chklim ) THEN
          WRITE (xern1,'(I8)') i
          WRITE (xern3,'(1PE15.6)') err
          CALL XERMSG('SLATEC','SNLS1','DERIVATIVE OF '//'FUNCTION '//&
            xern1//' MAY BE WRONG, ERR = '//xern3//&
            ' TOO CLOSE TO 0.',7,0)
        ENDIF
      ENDDO
      !
    ENDIF
  ENDIF
  IF ( iflag<0 ) GOTO 200
  !
  !        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
  !
  CALL QRFAC(M,N,Fjac,Ldfjac,.TRUE.,Ipvt,N,Wa1,Wa2,Wa3)
  !
  !        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
  !        QTF.
  !
  DO i = 1 , M
    Wa4(i) = Fvec(i)
  ENDDO
  DO j = 1 , N
    IF ( Fjac(j,j)/=zero ) THEN
      sum = zero
      DO i = j , M
        sum = sum + Fjac(i,j)*Wa4(i)
      ENDDO
      temp = -sum/Fjac(j,j)
      DO i = j , M
        Wa4(i) = Wa4(i) + Fjac(i,j)*temp
      ENDDO
    ENDIF
    Fjac(j,j) = Wa1(j)
    Qtf(j) = Wa4(j)
  ENDDO
ENDIF
!
!        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
IF ( iter==1 ) THEN
  IF ( Mode/=2 ) THEN
    DO j = 1 , N
      Diag(j) = Wa2(j)
      IF ( Wa2(j)==zero ) Diag(j) = one
    ENDDO
  ENDIF
  !
  !        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
  !        AND INITIALIZE THE STEP BOUND DELTA.
  !
  DO j = 1 , N
    Wa3(j) = Diag(j)*X(j)
  ENDDO
  xnorm = ENORM(N,Wa3)
  delta = Factor*xnorm
  IF ( delta==zero ) delta = Factor
ENDIF
!
!        COMPUTE THE NORM OF THE SCALED GRADIENT.
!
gnorm = zero
IF ( fnorm/=zero ) THEN
  DO j = 1 , N
    l = Ipvt(j)
    IF ( Wa2(l)/=zero ) THEN
      sum = zero
      DO i = 1 , j
        sum = sum + Fjac(i,j)*(Qtf(i)/fnorm)
      ENDDO
      gnorm = MAX(gnorm,ABS(sum/Wa2(l)))
    ENDIF
  ENDDO
ENDIF
!
!        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!
IF ( gnorm<=Gtol ) Info = 4
IF ( Info==0 ) THEN
  !
  !        RESCALE IF NECESSARY.
  !
  IF ( Mode/=2 ) THEN
    DO j = 1 , N
      Diag(j) = MAX(Diag(j),Wa2(j))
    ENDDO
  ENDIF
  DO
    !
    !        BEGINNING OF THE INNER LOOP.
    !
    !
    !           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
    !
    CALL LMPAR(N,Fjac,Ldfjac,Ipvt,Diag,Qtf,delta,par,Wa1,Wa2,Wa3,Wa4)
    !
    !           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
    !
    DO j = 1 , N
      Wa1(j) = -Wa1(j)
      Wa2(j) = X(j) + Wa1(j)
      Wa3(j) = Diag(j)*Wa1(j)
    ENDDO
    pnorm = ENORM(N,Wa3)
    !
    !           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
    !
    IF ( iter==1 ) delta = MIN(delta,pnorm)
    !
    !           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
    !
    iflag = 1
    CALL FCN(iflag,M,N,Wa2,Wa4,Fjac,ijunk)
    Nfev = Nfev + 1
    IF ( iflag<0 ) EXIT
    fnorm1 = ENORM(M,Wa4)
    !
    !           COMPUTE THE SCALED ACTUAL REDUCTION.
    !
    actred = -one
    IF ( p1*fnorm1<fnorm ) actred = one - (fnorm1/fnorm)**2
    !
    !           COMPUTE THE SCALED PREDICTED REDUCTION AND
    !           THE SCALED DIRECTIONAL DERIVATIVE.
    !
    DO j = 1 , N
      Wa3(j) = zero
      l = Ipvt(j)
      temp = Wa1(l)
      DO i = 1 , j
        Wa3(i) = Wa3(i) + Fjac(i,j)*temp
      ENDDO
    ENDDO
    temp1 = ENORM(N,Wa3)/fnorm
    temp2 = (SQRT(par)*pnorm)/fnorm
    prered = temp1**2 + temp2**2/p5
    dirder = -(temp1**2+temp2**2)
    !
    !           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
    !           REDUCTION.
    !
    ratio = zero
    IF ( prered/=zero ) ratio = actred/prered
    !
    !           UPDATE THE STEP BOUND.
    !
    IF ( ratio<=p25 ) THEN
      IF ( actred>=zero ) temp = p5
      IF ( actred<zero ) temp = p5*dirder/(dirder+p5*actred)
      IF ( p1*fnorm1>=fnorm.OR.temp<p1 ) temp = p1
      delta = temp*MIN(delta,pnorm/p1)
      par = par/temp
    ELSEIF ( par==zero.OR.ratio>=p75 ) THEN
      delta = pnorm/p5
      par = p5*par
    ENDIF
    !
    !           TEST FOR SUCCESSFUL ITERATION.
    !
    IF ( ratio>=p0001 ) THEN
      !
      !           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
      !
      DO j = 1 , N
        X(j) = Wa2(j)
        Wa2(j) = Diag(j)*X(j)
      ENDDO
      DO i = 1 , M
        Fvec(i) = Wa4(i)
      ENDDO
      xnorm = ENORM(N,Wa2)
      fnorm = fnorm1
      iter = iter + 1
    ENDIF
    !
    !           TESTS FOR CONVERGENCE.
    !
    IF ( ABS(actred)<=Ftol.AND.prered<=Ftol.AND.p5*ratio<=one ) Info = 1
    IF ( delta<=Xtol*xnorm ) Info = 2
    IF ( ABS(actred)<=Ftol.AND.prered<=Ftol.AND.p5*ratio<=one.AND.&
      Info==2 ) Info = 3
    IF ( Info/=0 ) EXIT
    !
    !           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
    !
    IF ( Nfev>=Maxfev ) Info = 5
    IF ( ABS(actred)<=epsmch.AND.prered<=epsmch.AND.p5*ratio<=one )&
      Info = 6
    IF ( delta<=epsmch*xnorm ) Info = 7
    IF ( gnorm<=epsmch ) Info = 8
    IF ( Info/=0 ) EXIT
    !
    !           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
    !
    !
    !        END OF THE OUTER LOOP.
    !
    IF ( ratio>=p0001 ) GOTO 100
  ENDDO
ENDIF
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
200  IF ( iflag<0 ) Info = iflag
iflag = 0
IF ( Nprint>0 ) CALL FCN(iflag,M,N,X,Fvec,Fjac,ijunk)
IF ( Info<0 ) CALL XERMSG('SLATEC','SNLS1',&
  'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.'&
  ,1,1)
IF ( Info==0 ) CALL XERMSG('SLATEC','SNLS1','INVALID INPUT PARAMETER.',2,&
  1)
IF ( Info==4 ) CALL XERMSG('SLATEC','SNLS1',&
  'THIRD CONVERGENCE CONDITION, CHECK RESULTS BEFORE ACCEPTING.'&
  ,1,1)
IF ( Info==5 ) CALL XERMSG('SLATEC','SNLS1',&
  'TOO MANY FUNCTION EVALUATIONS.',9,1)
IF ( Info>=6 ) CALL XERMSG('SLATEC','SNLS1',&
  'TOLERANCES TOO SMALL, NO FURTHER IMPROVEMENT POSSIBLE.'&
  ,3,1)
!
!     LAST CARD OF SUBROUTINE SNLS1.
!
END SUBROUTINE SNLS1
