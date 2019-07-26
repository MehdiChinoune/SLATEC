!** SNLS1E
PURE SUBROUTINE SNLS1E(FCN,Iopt,M,N,X,Fvec,Tol,Nprint,Info,Iw,Wa,Lwa)
  !> An easy-to-use code which minimizes the sum of the squares of M nonlinear functions
  !  in N variables by a modification of the Levenberg-Marquardt algorithm.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  K1B1A1, K1B1A2
  !***
  ! **Type:**      SINGLE PRECISION (SNLS1E-S, DNLS1E-D)
  !***
  ! **Keywords:**  EASY-TO-USE, LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
  !             NONLINEAR LEAST SQUARES
  !***
  ! **Author:**  Hiebert, K. L., (SNLA)
  !***
  ! **Description:**
  !
  ! 1. Purpose.
  !
  !       The purpose of SNLS1E is to minimize the sum of the squares of M
  !       nonlinear functions in N variables by a modification of the
  !       Levenberg-Marquardt algorithm.  This is done by using the more
  !       general least-squares solver SNLS1.  The user must provide a
  !       subroutine which calculates the functions.  The user has the
  !       option of how the Jacobian will be supplied.  The user can
  !       supply the full Jacobian, or the rows of the Jacobian (to avoid
  !       storing the full Jacobian), or let the code approximate the
  !       Jacobian by forward-differencing.  This code is the combination
  !       of the MINPACK codes (Argonne) LMDER1, LMDIF1, and LMSTR1.
  !
  !
  ! 2. Subroutine and Type Statements.
  !
  !       SUBROUTINE SNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
  !      *                  INFO,IW,WA,LWA)
  !       INTEGER IOPT,M,N,NPRINT,INFO,LWA
  !       INTEGER IW(N)
  !       REAL TOL
  !       REAL X(N),FVEC(M),WA(LWA)
  !       EXTERNAL FCN
  !
  !
  ! 3. Parameters.
  !
  !       Parameters designated as input parameters must be specified on
  !       entry to SNLS1E and are not changed on exit, while parameters
  !       designated as output parameters need not be specified on entry
  !       and are set to appropriate values on exit from SNLS1E.
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
  !         FJAC and LDFJAC may be ignored    , if IOPT=1.
  !         REAL FJAC(LDFJAC,N)               , if IOPT=2.
  !         REAL FJAC(N)                      , if IOPT=3.
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
  !         user wants to terminate execution of SNLS1E.  In this case,
  !         set IFLAG to a negative integer.
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
  !       M is a positive integer input variable set to the number of functions.
  !
  !       N is a positive integer input variable set to the number of variables.
  !         N must not exceed M.
  !
  !       X is an array of length N.  On input, X must contain an initial
  !         estimate of the solution vector.  On output, X contains the
  !         final estimate of the solution vector.
  !
  !       FVEC is an output array of length M which contains the functions
  !         evaluated at the output X.
  !
  !       TOL is a non-negative input variable.  Termination occurs when
  !         the algorithm estimates either that the relative error in the
  !         sum of squares is at most TOL or that the relative error
  !         between X and the solution is at most TOL.  Section 4 contains
  !         more details about TOL.
  !
  !       NPRINT is an integer input variable that enables controlled
  !         printing of iterates if it is positive.  In this case, FCN is
  !         called with IFLAG = 0 at the beginning of the first iteration
  !         and every NPRINT iterations thereafter and immediately prior
  !         to return, with X and FVEC available for printing. Appropriate
  !         print statements must be added to FCN (see example) and
  !         FVEC should not be altered.  If NPRINT is not positive, no
  !         special calls of FCN with IFLAG = 0 are made.
  !
  !       INFO is an integer output variable.  If the user has terminated
  !         execution, INFO is set to the (negative) value of IFLAG.  See
  !         description of FCN and JAC. Otherwise, INFO is set as follows.
  !
  !         INFO = 0  improper input parameters.
  !
  !         INFO = 1  algorithm estimates that the relative error in the
  !                   sum of squares is at most TOL.
  !
  !         INFO = 2  algorithm estimates that the relative error between
  !                   X and the solution is at most TOL.
  !
  !         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
  !
  !         INFO = 4  FVEC is orthogonal to the columns of the Jacobian to
  !                   machine precision.
  !
  !         INFO = 5  number of calls to FCN has reached 100*(N+1)
  !                   for IOPT=2 or 3 or 200*(N+1) for IOPT=1.
  !
  !         INFO = 6  TOL is too small.  No further reduction in the sum
  !                   of squares is possible.
  !
  !         INFO = 7  TOL is too small.  No further improvement in the
  !                   approximate solution X is possible.
  !
  !         Sections 4 and 5 contain more details about INFO.
  !
  !       IW is an INTEGER work array of length N.
  !
  !       WA is a work array of length LWA.
  !
  !       LWA is a positive integer input variable not less than
  !         N*(M+5)+M for IOPT=1 and 2 or N*(N+5)+M for IOPT=3.
  !
  !
  ! 4. Successful Completion.
  !
  !       The accuracy of SNLS1E is controlled by the convergence parame-
  !       ter TOL.  This parameter is used in tests which make three types
  !       of comparisons between the approximation X and a solution XSOL.
  !       SNLS1E terminates when any of the tests is satisfied.  If TOL is
  !       less than the machine precision (as defined by the function
  !       eps_sp), then SNLS1E only attempts to satisfy the test
  !       defined by the machine precision.  Further progress is not usu-
  !       ally possible.  Unless high precision solutions are required,
  !       the recommended value for TOL is the square root of the machine
  !       precision.
  !
  !       The tests assume that the functions are reasonably well behaved,
  !       and, if the Jacobian is supplied by the user, that the functions
  !       and the Jacobian are coded consistently.  If these conditions
  !       are not satisfied, then SNLS1E may incorrectly indicate conver-
  !       gence.  If the Jacobian is coded correctly or IOPT=1,
  !       then the validity of the answer can be checked, for example, by
  !       rerunning SNLS1E with tighter tolerances.
  !
  !       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
  !         of a vector Z, then this test attempts to guarantee that
  !
  !               ENORM(FVEC) <= (1+TOL)*ENORM(FVECS),
  !
  !         where FVECS denotes the functions evaluated at XSOL.  If this
  !         condition is satisfied with TOL = 10**(-K), then the final
  !         residual norm ENORM(FVEC) has K significant decimal digits and
  !         INFO is set to 1 (or to 3 if the second test is also satis-
  !         fied).
  !
  !       Second Convergence Test.  If D is a diagonal matrix (implicitly
  !         generated by SNLS1E) whose entries contain scale factors for
  !         the variables, then this test attempts to guarantee that
  !
  !               ENORM(D*(X-XSOL)) <=  TOL*ENORM(D*XSOL).
  !
  !         If this condition is satisfied with TOL = 10**(-K), then the
  !         larger components of D*X have K significant decimal digits and
  !         INFO is set to 2 (or to 3 if the first test is also satis-
  !         fied).  There is a danger that the smaller components of D*X
  !         may have large relative errors, but the choice of D is such
  !         that the accuracy of the components of X is usually related to
  !         their sensitivity.
  !
  !       Third Convergence Test.  This test is satisfied when FVEC is
  !         orthogonal to the columns of the Jacobian to machine preci-
  !         sion.  There is no clear relationship between this test and
  !         the accuracy of SNLS1E, and furthermore, the test is equally
  !         well satisfied at other critical points, namely maximizers and
  !         saddle points.  Therefore, termination caused by this test
  !         (INFO = 4) should be examined carefully.
  !
  !
  ! 5. Unsuccessful Completion.
  !
  !       Unsuccessful termination of SNLS1E can be due to improper input
  !       parameters, arithmetic interrupts, or an excessive number of
  !       function evaluations.
  !
  !       Improper Input Parameters.  INFO is set to 0 if IOPT < 1
  !         or IOPT > 3, or N <= 0, or M < N, or TOL < 0.E0,
  !         or for IOPT=1 or 2 LWA < N*(M+5)+M, or for IOPT=3
  !         LWA < N*(N+5)+M.
  !
  !       Arithmetic Interrupts.  If these interrupts occur in the FCN
  !         subroutine during an early stage of the computation, they may
  !         be caused by an unacceptable choice of X by SNLS1E.  In this
  !         case, it may be possible to remedy the situation by not evalu-
  !         ating the functions here, but instead setting the components
  !         of FVEC to numbers that exceed those in the initial FVEC.
  !
  !       Excessive Number of Function Evaluations.  If the number of
  !         calls to FCN reaches 100*(N+1) for IOPT=2 or 3 or 200*(N+1)
  !         for IOPT=1, then this indicates that the routine is converging
  !         very slowly as measured by the progress of FVEC, and INFO is
  !         set to 5.  In this case, it may be helpful to restart SNLS1E,
  !         thereby forcing it to disregard old (and possibly harmful)
  !         information.
  !
  !
  ! 6. Characteristics of the Algorithm.
  !
  !       SNLS1E is a modification of the Levenberg-Marquardt algorithm.
  !       Two of its main characteristics involve the proper use of
  !       implicitly scaled variables and an optimal choice for the cor-
  !       rection.  The use of implicitly scaled variables achieves scale
  !       invariance of SNLS1E and limits the size of the correction in
  !       any direction where the functions are changing rapidly.  The
  !       optimal choice of the correction guarantees (under reasonable
  !       conditions) global convergence from starting points far from the
  !       solution and a fast rate of convergence for problems with small
  !       residuals.
  !
  !       Timing.  The time required by SNLS1E to solve a given problem
  !         depends on M and N, the behavior of the functions, the accu-
  !         racy requested, and the starting point.  The number of arith-
  !         metic operations needed by SNLS1E is about N**3 to process
  !         each evaluation of the functions (call to FCN) and to process
  !         each evaluation of the Jacobian SNLS1E takes M*N**2 for IOPT=2
  !         (one call to JAC), M*N**2 for IOPT=1 (N calls to FCN) and
  !         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
  !         can be evaluated quickly, the timing of SNLS1E will be
  !         strongly influenced by the time spent in FCN.
  !
  !       Storage.  SNLS1E requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
  !         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
  !         locations and N integer storage locations, in addition to
  !         the storage required by the program.  There are no internally
  !         declared storage arrays.
  !
  !- Long Description:
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
  ! C     Driver for SNLS1E example.
  ! C
  !       INTEGER I,IOPT,M,N,NPRINT,JNFO,LWA,NWRITE
  !       INTEGER IW(3)
  !       REAL TOL,FNORM
  !       REAL X(3),FVEC(15),WA(75)
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
  !       LWA = 75
  !       NPRINT = 0
  ! C
  ! C     Set TOL to the square root of the machine precision.
  ! C     Unless high precision solutions are required,
  ! C     this is the recommended setting.
  ! C
  !       TOL = SQRT(eps_sp)
  ! C
  !       CALL SNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
  !      *            INFO,IW,WA,LWA)
  !       FNORM = ENORM(M,FVEC)
  !       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
  !       STOP
  !  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
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
  !       IF(IFLAG /= 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF(I > 8) TMP3 = TMP2
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
  !       IF(IFLAG /= 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       IF(IFLAG/=1) GO TO 20
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF(I > 8) TMP3 = TMP2
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
  !          IF(I > 8) TMP3 = TMP2
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
  !       IF(IFLAG /= 0) GO TO 5
  ! C
  ! C     Insert print statements here when NPRINT is positive.
  ! C
  !       RETURN
  !     5 CONTINUE
  !       IF( IFLAG/=1) GO TO 20
  !       DO 10 I = 1, M
  !          TMP1 = I
  !          TMP2 = 16 - I
  !          TMP3 = TMP1
  !          IF(I > 8) TMP3 = TMP2
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
  !          IF(I > 8) TMP3 = TMP2
  !          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
  !          FJAC(1) = -1.E0
  !          FJAC(2) = TMP1*TMP2/TMP4
  !          FJAC(3) = TMP1*TMP3/TMP4
  !       RETURN
  !       END
  !
  !***
  ! **References:**  Jorge J. More, The Levenberg-Marquardt algorithm:
  !                 implementation and theory.  In Numerical Analysis
  !                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
  !                 Watson, Editor), Lecture Notes in Mathematics 630,
  !                 Springer-Verlag, 1978.
  !***
  ! **Routines called:**  SNLS1, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTERFACE
    PURE SUBROUTINE FCN(Iflag,M,N,X,Fvec,Fjac,Ldfjac)
      IMPORT SP
      INTEGER, INTENT(IN) :: Ldfjac, M, N, Iflag
      REAL(SP), INTENT(IN) :: X(N)
      REAL(SP), INTENT(INOUT) :: Fvec(M)
      REAL(SP), INTENT(OUT) :: Fjac(:,:)
    END SUBROUTINE FCN
  END INTERFACE
  INTEGER, INTENT(IN) :: M, N, Nprint, Lwa, Iopt
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(OUT) :: Iw(N)
  REAL(SP), INTENT(IN) :: Tol
  REAL(SP), INTENT(INOUT) :: X(N)
  REAL(SP), INTENT(OUT) :: Fvec(M), Wa(Lwa)
  !
  INTEGER :: indexx, maxfev, mode, nfev, njev
  REAL(SP) :: ftol, gtol, xtol, epsfcn
  REAL(SP), PARAMETER :: factor = 100._SP
  !* FIRST EXECUTABLE STATEMENT  SNLS1E
  Info = 0
  !
  !     CHECK THE INPUT PARAMETERS FOR ERRORS.
  !
  IF( Iopt>=1 .AND. Iopt<=3 .AND. N>0 .AND. M>=N .AND. Tol>=0._SP .AND. Lwa>=N*(N+5)+M ) THEN
    IF( Iopt>=3 .OR. Lwa>=N*(M+5)+M ) THEN
      !  CALL SNLS1.
      maxfev = 100*(N+1)
      IF( Iopt==1 ) maxfev = 2*maxfev
      ftol = Tol
      xtol = Tol
      gtol = 0._SP
      epsfcn = 0._SP
      mode = 1
      indexx = 5*N + M
      CALL SNLS1(FCN,Iopt,M,N,X,Fvec,Wa(indexx+1),M,ftol,xtol,gtol,maxfev,&
        epsfcn,Wa(1),mode,factor,Nprint,Info,nfev,njev,Iw,Wa(N+1),&
        Wa(2*N+1),Wa(3*N+1),Wa(4*N+1),Wa(5*N+1))
      IF( Info==8 ) Info = 4
    END IF
  END IF
  IF( Info==0 ) ERROR STOP 'SNLS1E : INVALID INPUT PARAMETER.'
  !
  !     LAST CARD OF SUBROUTINE SNLS1E.
  !
END SUBROUTINE SNLS1E